# Joint Calling Workflow
# Copyright (c) 2019 Genome Research Limited
# Licensed under GPLv3, or later

# Adapted from exome joint calling pipeline from the Broad Institute
# Copyright (c) 2019 Broad Institute
# Licensed under the BSD 3-Clause License

workflow JointCalling {
  # Inputs
  # n.b., Use camelCase for clarity; snake_case is used for derivatives
  File referenceFASTA     # Reference FASTA file
  File referenceIndex     # Reference index
  File referenceDict      # Reference dictionary
  File sampleNameMap      # Sample name mapping
  File unpaddedIntervals  # List of unpadded intervals
  File dbsnpVCF           # dbSNP VCF
  Int? vcfCount           # Optional: Scatter count

  Array[Array[String]] sample_name_map_lines = read_tsv(sampleNameMap)
  Int gvcf_count = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval
  # list. We allow overriding the behaviour by specifying the desired
  # number of VCFs to scatter over for testing/special requests.
  # n.b., WGS runs get 30x more scattering than exome and exome scatter
  # count per sample is 0.05 (modulo 10 <= scatter_count <= 1000)
  # TODO Bound scatter_count by 1000
  # TODO If vcfCount is specified, then presumably it must be <= gvcf_count
  Int unbounded_scatter_count = select_first([vcfCount, round(0.05 * gvcf_count)])
  Int scatter_count = if unbounded_scatter_count > 10 then unbounded_scatter_count else 10

  call SplitIntervalList {
    input:
      referenceFASTA = referenceFASTA,
      referenceIndex = referenceIndex,
      referenceDict  = referenceDict,
      intervalList   = unpaddedIntervals,
      scatter_count  = scatter_count
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it is the
    # optimal value for the amount of memory allocated within the task;
    # please do not change it without consulting the GATK engine team!
    call ImportGVCFs {
      input:
        sampleNameMap      = sampleNameMap,
        interval           = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        batch_size         = 50
    }

    call GenotypeGVCFs {
      input:
        referenceFASTA      = referenceFASTA,
        referenceIndex      = referenceIndex,
        referenceDict       = referenceDict,
        dbsnpVCF            = dbsnpVCF,
        workspace_tar       = ImportGVCFs.output_genomicsdb,
        interval            = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz"
    }
  }

  # TODO We ought to have a gather stage here, surely? Broad's original
  # task uses GatherVcfsCloud in their VQSR step; we'll have to use
  # vanilla GatherVcfs...

  output {
    # TODO Gathered output from scatter
  }
}

task SplitIntervalList {
  # The raw intervals file given as input to the workflow might contain
  # a TON of entries, to the point that scattering over all of them
  # would kill Cromwell. This task scans the file for contiguous
  # intervals which can be processed together to reduce the scatter
  # width. For example:
  #
  #   chr1:1-195878       \      chr1:1-391754
  #   chr1:195879-391754   \___  chr2:1-161787
  #   chr2:1-161787        /     chr2:323574-323584
  #   chr2:323574-323584  /

  String intervalList
  Int    scatter_count
  File   referenceFASTA
  File   referenceIndex
  File   referenceDict

  command <<<
    /gatk/gatk SplitIntervals \
      -L "${intervalList}" -O scatterDir -scatter ${scatter_count} -R "${referenceFASTA}" \
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
  >>>

  runtime {
    lsf_memory:  3072

    # TODO We use Laura Gauthier's GATK fork for joint calling
    # (4.0.11.0-22-gae8e9f0-SNAPSHOT), which we've pressed into a
    # Singularity 2.5.2 image. Move this to production GATK in a SIF
    # (Singularity 3.0) container once the necessary conditions are met.
    singularity: "/software/hgi/containers/gatk-jointcalling.simg"
    # singularity: "/software/hgi/containers/gatk-4.1.0.0.sif"
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

task ImportGVCFs {
  File   sampleNameMap
  File   interval
  String workspace_dir_name
  Int    batch_size

  command <<<
    set -euo pipefail

    rm -rf "${workspace_dir_name}"

    # We've seen some GenomicsDB performance regressions related to
    # intervals, so we're going to only supply a single interval.
    # There's no data in between since we didn't run HaplotypeCaller
    # over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GB
    # lower than the total allocation because this tool uses a
    # significant amount of non-heap memory for native libraries. Also,
    # testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads.
    # TODO Set -Xms option based on lsf_memory, rather than hardcoded
    /gatk/gatk --java-options -Xms4g \
      GenomicsDBImport \
      --genomicsdb-workspace-path "${workspace_dir_name}" \
      --batch-size ${batch_size} \
      -L "${interval}" \
      --sample-name-map "${sampleNameMap}" \
      --reader-threads 1 \
      -ip 500

    tar -cf "${workspace_dir_name}.tar" "${workspace_dir_name}"
  >>>

  runtime {
    lsf_memory:  7168
    lsf_cores:   2

    # TODO We use Laura Gauthier's GATK fork for joint calling
    # (4.0.11.0-22-gae8e9f0-SNAPSHOT), which we've pressed into a
    # Singularity 2.5.2 image. Move this to production GATK in a SIF
    # (Singularity 3.0) container once the necessary conditions are met.
    singularity: "/software/hgi/containers/gatk-jointcalling.simg"
    # singularity: "/software/hgi/containers/gatk-4.1.0.0.sif"
  }

  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File   referenceFASTA
  File   referenceIndex
  File   referenceDict
  File   dbsnpVCF
  File   workspace_tar
  String interval
  String output_vcf_filename

  command <<<
    set -euo pipefail

    tar -xf "${workspace_tar}"
    WORKSPACE="$(basename "${workspace_tar}" .tar)"

    /gatk/gatk SpanIntervals \
      -L "${interval}" -R "${referenceFASTA}" -O spanning.interval_list

    # TODO Set -Xms option based on lsf_memory, rather than hardcoded
    # (see note in ImportGVCFs task, above, for justification)
    /gatk/gatk --java-options -Xms5g \
      GenotypeGVCFs \
      -R "${referenceFASTA}" \
      -O "${output_vcf_filename}" \
      -D "${dbsnpVCF}" \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      --use-new-qual-calculator \
      -V "gendb://$WORKSPACE" \
      -L spanning.interval_list
  >>>

  runtime {
    lsf_memory:  7168
    lsf_cores:   2

    # TODO We use Laura Gauthier's GATK fork for joint calling
    # (4.0.11.0-22-gae8e9f0-SNAPSHOT), which we've pressed into a
    # Singularity 2.5.2 image. Move this to production GATK in a SIF
    # (Singularity 3.0) container once the necessary conditions are met.
    singularity: "/software/hgi/containers/gatk-jointcalling.simg"
    # singularity: "/software/hgi/containers/gatk-4.1.0.0.sif"
  }

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}
