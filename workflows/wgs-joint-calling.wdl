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


  Int num_of_original_intervals = length(read_lines(unpaddedIntervals))
  Int num_gvcfs = length(read_lines(sampleNameMap))

  # Make a 2.5:1 interval number to samples in callset ratio interval list
  Int possible_merge_count = floor(num_of_original_intervals / num_gvcfs / 2.5)
  Int merge_count = if possible_merge_count > 1 then possible_merge_count else 1

  call DynamicallyCombineIntervals {
    input:
      intervals = unpaddedIntervals,
      merge_count = merge_count
  }

  Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)

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

task ImportGVCFs {
  File sampleNameMap
  String interval
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
    /gatk/gatk --java-options "-Xms15g -Xmx15g" \
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
    lsf_memory: 20720
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
  String interval
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
    /gatk/gatk --java-options "-Xms15g -Xmx15g" \
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
    lsf_memory: 20720
    lsf_cores:  2
    lsf_queue: "long"

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

task DynamicallyCombineIntervals {
  File intervals
  Int merge_count

  command {
    python << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ${merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("${intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
  }

  runtime {
    memory: "3 GB"
    preemptible: 5
    docker: "python:2.7"
  }

  output {
    File output_intervals = "out.intervals"
  }
}
