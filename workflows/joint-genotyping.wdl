# Copyright (c) 2019 Genome Research Limited
# Licensed under GPLv3, or later

# Adapted from exome joint calling pipeline from the Broad Institute
# Copyright (c) 2019 Broad Institute
# Licensed under the BSD 3-Clause License

workflow JointGenotypingForExomes {
  File unpadded_intervals_file

  String callset_name
  File sample_name_map

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbsnp_vcf
  File dbsnp_vcf_index

  Int small_disk
  Int medium_disk

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000
  Int? vcf_count
  Int unboundedScatterCount = select_first([vcf_count, round(0.05 * num_gvcfs)])
  Int scatterCount = if unboundedScatterCount > 10 then unboundedScatterCount else 10 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?
  call SplitIntervalList {
    input:
      intervalList = unpadded_intervals_file,
      scatterCount = scatterCount,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        disk_size = medium_disk,
        batch_size = 50
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        disk_size = medium_disk
    }
  }

  output {
    # TODO
  }
}

## The raw intervals file given as input to the workflow might contain a TON of entries, to
## the point that scattering over all of them would kill Cromwell. This task scans the file
## for contiguous intervals which can be processed together to reduce the scatter width.
##
## Input intervals are expected to look like:
##     chr1:1-195878
##     chr1:195879-391754
##     chr2:1-161787
##     chr2:323574-323584
##
## For that specific input, the output would look like:
##     chr1:1-391754
##     chr2:1-161787
##     chr2:323574-323584
#

task SplitIntervalList {
  String intervalList
  Int disk_size
  Int scatterCount
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  command <<<
    /usr/gitc/gatk SplitIntervals \
      -L ${intervalList} -O  scatterDir -scatter ${scatterCount} -R ${ref_fasta} \
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
   >>>

  runtime {
    memory: "3 GB"
    preemptible: 5
    disks: "local-disk " + disk_size + " HDD"
    docker: "ldgauthier/gatk_exome_joint_calling"
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

task ImportGVCFs {
  File sample_name_map
  File interval
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String workspace_dir_name

  Int disk_size
  Int batch_size

  command <<<
    set -euo pipefail

    rm -rf ${workspace_dir_name}

    #We've seen some GenomicsDB performance regressions related to intervals, so we're going to only supply a single interval
    #There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute
    #/usr/gitc/gatk SpanIntervals -L ${interval} -R ${ref_fasta} -O spanning.interval_list
    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    /usr/gitc/gatk --java-options -Xms4g \
      GenomicsDBImport \
      --genomicsdb-workspace-path ${workspace_dir_name} \
      --batch-size ${batch_size} \
      -L ${interval} \
      --sample-name-map ${sample_name_map} \
      --reader-threads 1 \
      -ip 500

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}
  >>>

  runtime {
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
    docker: "ldgauthier/gatk_exome_joint_calling"
  }

  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String dbsnp_vcf

  Int disk_size

  command <<<
    set -euo pipefail

    tar -xf ${workspace_tar}
    WORKSPACE=$(basename ${workspace_tar} .tar)

    /usr/gitc/gatk SpanIntervals -L ${interval} -R ${ref_fasta} -O spanning.interval_list

    /usr/gitc/gatk --java-options -Xms5g \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbsnp_vcf} \
     -G StandardAnnotation -G AS_StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://$WORKSPACE \
     -L spanning.interval_list
  >>>

  runtime {
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
    docker: "ldgauthier/gatk_exome_joint_calling"
  }

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}
