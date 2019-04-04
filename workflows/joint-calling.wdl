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
  File dbsnpVCFIndex      # dbSNP VCF index
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
        sampleNameMap = sampleNameMap,
        interval      = unpadded_intervals[idx],
        batch_size    = 50
    }

    call GenotypeGVCFs {
      input:
        referenceFASTA      = referenceFASTA,
        referenceIndex      = referenceIndex,
        referenceDict       = referenceDict,
        dbsnpVCF            = dbsnpVCF,
        dbsnpVCFIndex       = dbsnpVCFIndex,
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
    /gatk/gatk --java-options "-Xms3g -Xmx3g" \
      SplitIntervals \
      -L "${intervalList}" -O scatterDir -scatter ${scatter_count} -R "${referenceFASTA}" \
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
  >>>

  runtime {
    lsf_memory: 3072

    # TODO Move this to SIF (Singularity 3.0) container when ready
    singularity: "/software/hgi/containers/gatk-4.1.0.0.simg"
    # singularity: "/software/hgi/containers/gatk-4.1.0.0.sif"
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

task ImportGVCFs {
  File sampleNameMap
  File interval
  Int  batch_size

  command <<<
    set -eu

    # TODO Make TMPDIR parametrisable; this would best be done via a
    # runtime attribute, but tasks apparently don't have any visibility
    # of these in Cromwell.
    # FIXME Despite the documentation, GATK doesn't like it if the
    # temporary directory already exists (even when it's empty); the
    # below workaround is subject to race conditions
    declare WORKSPACE="$(TMPDIR="/tmp" mktemp -du)"
    trap 'rm -rf "$WORKSPACE"' EXIT

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
    /gatk/gatk --java-options "-Xms25g -Xmx25g" \
      GenomicsDBImport \
      --genomicsdb-workspace-path "$WORKSPACE" \
      --batch-size ${batch_size} \
      -L "${interval}" \
      --merge-input-intervals \
      --sample-name-map "${sampleNameMap}" \
      --reader-threads 1 \
      -ip 500

    tar cf "genomicsdb.tar" -C "$WORKSPACE" .
  >>>

  runtime {
    lsf_memory: 30720
    lsf_cores:  2

    # TODO Temporary space requirement should be a function of,
    # presumably, the interval size; set to 20GiB for now...
    lsf_resources: "rusage[tmp=20480]"

    # TODO Move this to SIF (Singularity 3.0) container when ready
    singularity: "/software/hgi/containers/gatk-4.1.0.0.simg"
    # singularity: "/software/hgi/containers/gatk-4.1.0.0.sif"
  }

  output {
    File output_genomicsdb = "genomicsdb.tar"
  }
}

task GenotypeGVCFs {
  File   referenceFASTA
  File   referenceIndex
  File   referenceDict
  File   dbsnpVCF
  File   dbsnpVCFIndex
  File   workspace_tar
  File   interval
  String output_vcf_filename

  command <<<
    set -eu

    # TODO Make TMPDIR parametrisable; this would best be done via a
    # runtime attribute, but tasks apparently don't have any visibility
    # of these in Cromwell.
    declare WORKSPACE="$(TMPDIR="/tmp" mktemp -d)"
    trap 'rm -rf "$WORKSPACE"' EXIT
    tar xf "${workspace_tar}" -C "$WORKSPACE"

    # TODO Set -Xms option based on lsf_memory, rather than hardcoded
    # (see note in ImportGVCFs task, above, for justification)
    /gatk/gatk --java-options "-Xms5g -Xmx5g" \
      GenotypeGVCFs \
      -R "${referenceFASTA}" \
      -O "${output_vcf_filename}" \
      -D "${dbsnpVCF}" \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      --use-new-qual-calculator \
      -V "gendb://$WORKSPACE" \
      -L "${interval}"
  >>>

  runtime {
    lsf_memory: 7168
    lsf_cores:  2

    # TODO Temporary space requirement should be a function of,
    # presumably, the interval size; set to 20GiB for now...
    lsf_resources: "rusage[tmp=20480]"

    # TODO Move this to SIF (Singularity 3.0) container when ready
    singularity: "/software/hgi/containers/gatk-4.1.0.0.simg"
    # singularity: "/software/hgi/containers/gatk-4.1.0.0.sif"
  }

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}
