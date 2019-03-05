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
  Int large_disk
  Int huge_disk

  Array[String] snp_recalibration_tranche_values
  Array[String] snp_recalibration_annotation_values
  Array[String] indel_recalibration_tranche_values
  Array[String] indel_recalibration_annotation_values

  File haplotype_database

  File eval_interval_list
  File hapmap_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf
  File one_thousand_genomes_resource_vcf_index
  File mills_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf = dbsnp_vcf
  File dbsnp_resource_vcf_index = dbsnp_vcf_index

  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69
  Float snp_filter_level
  Float indel_filter_level
  Int SNP_VQSR_downsampleFactor

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

    call HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        vcf_index = GenotypeGVCFs.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size = medium_disk
    }
  }

  call GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs_fofn = write_lines(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf),
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size = medium_disk
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      disk_size = small_disk
  }

  if (num_gvcfs > 500000) {
    call SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        downsampleFactor = SNP_VQSR_downsampleFactor,
        model_report_filename = callset_name + ".snps.model.report",
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        disk_size = small_disk
    }


    scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
      call SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
          sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
          recalibration_filename = callset_name + ".snps." + idx + ".recal",
          tranches_filename = callset_name + ".snps." + idx + ".tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          disk_size = small_disk
        }
    }

    call GatherTranches as SNPGatherTranches {
      input:
        input_fofn = write_lines(SNPsVariantRecalibratorScattered.tranches),
        output_filename = callset_name + ".snps.gathered.tranches",
        disk_size = small_disk
    }
  }

  if (num_gvcfs <= 10000){
    call SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        disk_size = small_disk
    }
  }

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  # We allow overriding this default behavior for testing / special requests.
  Boolean? gather_vcfs
  Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    call ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
        snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
        snps_tranches = select_first([SNPGatherTranches.tranches, SNPsVariantRecalibratorClassic.tranches]),
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        disk_size = medium_disk
    }

    # For large callsets we need to collect metrics from the shards and gather them later.
    if (!is_small_callset) {
      call CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size = medium_disk
      }
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
  if (is_small_callset) {
    call GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs_fofn = write_lines(ApplyRecalibration.recalibrated_vcf),
        output_vcf_name = callset_name + ".vcf.gz",
        disk_size = huge_disk
    }

    call CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size = large_disk
    }
  }

  if (!is_small_callset) {
    # For large callsets we still need to gather the sharded metrics.
    call GatherMetrics {
      input:
        input_details_fofn = write_lines(select_all(CollectMetricsSharded.detail_metrics_file)),
        input_summaries_fofn = write_lines(select_all(CollectMetricsSharded.summary_metrics_file)),
        output_prefix = callset_name,
        disk_size = medium_disk
    }

    # We also need to determine which of the sharded VCFs contain fingerprinting sites,
    # so the final call to CrossCheckFingerprints doesn't have to pointlessly download
    # thousands of files.
    call GetFingerprintSiteVcfs {
      input:
        vcf_names = ApplyRecalibration.recalibrated_vcf,
        intervals = SplitIntervalList.output_intervals,
        haplotype_database = haplotype_database
    }
  }

  Array[File] vcfs_for_crosscheck = if is_small_callset
    then select_all([FinalGatherVcf.output_vcf])
    else select_first([GetFingerprintSiteVcfs.fingerprint_site_vcfs])

  scatter (row in sample_name_map_lines) {
    # GenomicDbImport could rename samples, creating mismatches between
    # the input GVCFs and the generated VCFs. We need to correct the GVCFs
    # to match before feeding them into CrosscheckFingerprints.
    call RenameGvcfSample {
      input:
        gvcf = row[1],
        gvcf_index = row[1] + ".tbi",
        new_sample_name = row[0]
    }
  }

  call CrossCheckFingerprint {
    input:
      gvcf_paths = RenameGvcfSample.renamed_gvcf,
      vcf_paths = vcfs_for_crosscheck,
      haplotype_database = haplotype_database,
      output_base_name = callset_name
  }

  output {
    # Outputs from the small callset path through the wdl.
    FinalGatherVcf.output_vcf
    FinalGatherVcf.output_vcf_index
    CollectMetricsOnFullVcf.detail_metrics_file
    CollectMetricsOnFullVcf.summary_metrics_file

    # Outputs from the large callset path through the wdl.
    # (note that we do not list ApplyRecalibration here because it is run in both paths)
    GatherMetrics.detail_metrics_file
    GatherMetrics.summary_metrics_file

    # Output the interval list generated/used by this run workflow.
    SplitIntervalList.output_intervals

    # Output the metrics from crosschecking fingerprints.
    CrossCheckFingerprint.output_check
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
task DynamicallyCombineIntervals {
  File intervals
  Int target_interval_count

  String combined_intervals_output = "out.intervals"

  command <<<
    set -euo pipefail

    python <<CODE

    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    # List of (chr, start, end) for the non-contiguous intervals
    # described in the unpadded intervals file.
    distinct_intervals = []

    # chr, start, and end for the current contiguous interval.
    current_chr, current_start, current_end = None, 0, 0

    # First build a list of all the intervals in the unpadded file,
    # collapsing contiguous intervals together in the process.
    with open("${intervals}") as intervals:
        for interval_line in intervals.readlines():
            chr, start, end = parse_interval(interval_line)

            # Initialization.
            if not current_chr:
                current_chr, current_start, current_end = chr, start, end

            # Keep the chain going if adjacent.
            elif chr == current_chr and start == current_end + 1:
                current_end = end

            # Not adjacent, start a new chain.
            else:
                distinct_intervals.append((current_chr, current_start, current_end))
                current_chr, current_start, current_end = chr, start, end

        # Last-seen interval, nothing more to chain.
        distinct_intervals.append((current_chr, current_start, current_end))

    # If we haven't reached the target number of output intervals, split the largest
    # intervals in half until we reach the target count. There's a possibility that we'll
    # have MORE than the target number of intervals before reaching this loop. In that
    # case there's nothing we can do, since the input intervals had more non-contiguous
    # regions than our target number.
    while len(distinct_intervals) < ${target_interval_count}:
        sorted_intervals = sorted(distinct_intervals, key=lambda i: i[2] - i[1])
        chr, start, end = sorted_intervals.pop()

        interval_size = end - start + 1

        # Pathological case: all intervals of size 1, and not enough of them to reach
        # the target count.
        if interval_size <= 1:
            break

        l_start = start
        l_end = start + (interval_size / 2)
        r_start = l_end + 1
        r_end = end

        distinct_intervals = sorted_intervals + [(chr, l_start, l_end), (chr, r_start, r_end)]

    # Reassemble the final intervals list in the correct order, so downstream tools don't break.
    # This is hard-coded for hg38, and it expects chromosome names of form chr<N> for <N> in {1-22, X, Y, M}.
    def compare_intervals(l, r):
        l_chr, l_start, l_end = l
        r_chr, r_start, r_end = r

        chr_values = dict(("chr%d"%(i + 1), i + 1) for i in range(22))
        chr_values["chrX"] = 23
        chr_values["chrY"] = 24
        chr_values["chrM"] = 25

        chr_cmp = cmp(chr_values[l_chr], chr_values[r_chr])
        if chr_cmp == 0:
            return cmp(l_start, r_start)
        else:
            return chr_cmp

    distinct_intervals.sort(compare_intervals)

    lines_to_write = [chr + ":" + str(start) + "-" + str(end) for chr, start, end in distinct_intervals]
    with open("${combined_intervals_output}", "w") as out:
        out.writelines("\n".join(lines_to_write))
    CODE
  >>>

  runtime {
    memory: "3 GB"
    preemptible: 5
    docker: "python:2.7"
  }

  output {
    File output_intervals = combined_intervals_output
  }
}

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
    docker: "ldgauthier/exome_joint_calling@sha256:b234d6b1af4b0a5de5dda28c8516209986a95c579b5ebf75b0943e7850a226a5"
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
    docker: "ldgauthier/exome_joint_calling@sha256:b234d6b1af4b0a5de5dda28c8516209986a95c579b5ebf75b0943e7850a226a5"
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
    docker: "ldgauthier/exome_joint_calling@sha256:b234d6b1af4b0a5de5dda28c8516209986a95c579b5ebf75b0943e7850a226a5"
  }

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename

  Int disk_size

  command <<<
    set -euo pipefail

    /usr/gitc/gatk --java-options -Xms3g \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    java -Xms3g -jar /usr/gitc/picard.jar \
      MakeSitesOnlyVcf \
      INPUT=${variant_filtered_vcf_filename} \
      OUTPUT=${sites_only_vcf_filename}
  >>>

  runtime {
    memory: "3.5 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  }
}

task IndelsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File dbsnp_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command <<<
    set -euo pipefail

    /usr/gitc/gatk --java-options -Xms24g \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      --use-allele-specific-annotations \
      -mode INDEL \
      --max-gaussians 4 \
      -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "26 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task SNPsVariantRecalibratorCreateModel {
  String recalibration_filename
  String tranches_filename
  Int downsampleFactor
  String model_report_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command <<<
    set -euo pipefail

    /usr/gitc/gatk --java-options -Xms100g \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      --use-allele-specific-annotations \
      -mode SNP \
      --sample-every-Nth-variant ${downsampleFactor} \
      --output-model ${model_report_filename} \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "104 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File model_report = "${model_report_filename}"
  }
}

task SNPsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename
  File? model_report

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index

  Int disk_size

  command <<<
    set -euo pipefail

    /usr/gitc/gatk --java-options -Xms3g \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      --use-allele-specific-annotations \
      -mode SNP \
       ${"--input-model " + model_report + " --output-tranches-for-scatter "} \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "3.5 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task GatherTranches {
  File input_fofn
  String output_filename

  Int disk_size

  command <<<
    set -euo pipefail

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf tranches
    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat ${input_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I tranches/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    cat ${input_fofn} | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

    /usr/gitc/gatk --java-options -Xms6g \
      GatherTranches \
      --input inputs.list \
      --output ${output_filename}
  >>>

  runtime {
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File tranches = "${output_filename}"
  }
}

task ApplyRecalibration {
  String recalibrated_vcf_filename
  File input_vcf
  File input_vcf_index
  File indels_recalibration
  File indels_recalibration_index
  File indels_tranches
  File snps_recalibration
  File snps_recalibration_index
  File snps_tranches

  Float indel_filter_level
  Float snp_filter_level

  Int disk_size

  command <<<
    set -euo pipefail

    /usr/gitc/gatk --java-options -Xms5g \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --use-allele-specific-annotations \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    /usr/gitc/gatk --java-options -Xms5g \
      ApplyVQSR \
      -O ${recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snps_recalibration} \
      --use-allele-specific-annotations \
      --tranches-file ${snps_tranches} \
      --truth-sensitivity-filter-level ${snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  >>>

  runtime {
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}

task GatherVcfs {
  File input_vcfs_fofn
  String output_vcf_name

  Int disk_size

  command <<<
    set -euo pipefail

    # Now using NIO to localize the vcfs but the input file must have a ".list" extension
    mv ${input_vcfs_fofn} inputs.list

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    /usr/gitc/gatk --java-options -Xms6g \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input inputs.list \
      --output ${output_vcf_name}

    tabix ${output_vcf_name}
  >>>

  runtime {
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {
  File input_vcf
  File input_vcf_index

  String metrics_filename_prefix
  File dbsnp_vcf
  File dbsnp_vcf_index
  File interval_list
  File ref_dict

  Int disk_size

  command <<<
    set -euo pipefail

    java -Xms6g -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      DBSNP=${dbsnp_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
      OUTPUT=${metrics_filename_prefix} \
      THREAD_COUNT=8 \
      TARGET_INTERVALS=${interval_list}
  >>>

  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }

  runtime {
    memory: "7 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
}

task GatherMetrics {
  File input_details_fofn
  File input_summaries_fofn

  String output_prefix

  Int disk_size

  command <<<
    set -euo pipefail

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat ${input_details_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat ${input_summaries_fofn} | /root/google-cloud-sdk/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=$(cat ${input_details_fofn} | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("I=metrics/%s ", $1)}')

    java -Xms2g -jar /usr/gitc/picard.jar \
      AccumulateVariantCallingMetrics \
      $INPUT \
      O=${output_prefix}
  >>>

  runtime {
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}

task GetFingerprintSiteVcfs {
  Array[String] vcf_names

  File intervals
  File haplotype_database

  String intervals_bed = "intervals.bed"
  String haplotype_db_bed = "hapmap.bed"

  command <<<
    set -euo pipefail

    python <<CODE
    import csv

    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    with open("${intervals}") as intervals:
      with open("${intervals_bed}", "w") as bed:
        n = 0
        writer = csv.writer(bed, delimiter='\t')
        for line in intervals.readlines():
          chr, start, end = parse_interval(line)
          writer.writerow([chr, start - 1, end, n])
          n += 1

    with open("${haplotype_database}") as db:
      with open("${haplotype_db_bed}", "w") as bed:
        writer = csv.writer(bed, delimiter='\t')
        for line in db.readlines():
          if line.startswith('@') or line.startswith('#'):
            continue
          chr, loc = tuple(line.split()[:2])
          writer.writerow([chr, int(loc) - 1, loc])
    CODE

    /app/bedtools intersect -u -wa -a ${intervals_bed} -b ${haplotype_db_bed} | cut -f4 > vcf_indexes

    # Quoth StackExchange: https://unix.stackexchange.com/a/418272
    awk 'NR==FNR{ pos[$1]; next }FNR in pos' vcf_indexes ${write_lines(vcf_names)}
  >>>

  runtime {
    memory: "3 GB"
    preemtible: 5
    docker: "biocontainers/bedtools"
  }

  output {
    Array[String] fingerprint_site_vcfs = read_lines(stdout())
  }
}

task RenameGvcfSample {
  File gvcf
  File gvcf_index
  String new_sample_name

  Int disk_size = round(2 * (size(gvcf, "GB") + size(gvcf_index, "GB")) + 1)
  String renamed = "renamed_" + basename(gvcf)

  command <<<
    set -euo pipefail

    java -Xms2g -jar /usr/gitc/picard.jar \
      RenameSampleInVcf \
      I=${gvcf} \
      O=${renamed} \
      NEW_SAMPLE_NAME='${new_sample_name}' \
      CREATE_INDEX=true
  >>>

  runtime {
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }

  output {
    File renamed_gvcf = renamed
    File renamed_gvcf_index = renamed + ".tbi"
  }
}

task CrossCheckFingerprint {
  Array[String] gvcf_paths
  Array[String] vcf_paths
  File haplotype_database
  String output_base_name

  Int num_gvcfs = length(gvcf_paths)
  Int cpu = if num_gvcfs < 32 then num_gvcfs else 32
  # Compute memory to use based on the CPU count, following the pattern of
  # 3.75GB / cpu used by GCP's pricing: https://cloud.google.com/compute/pricing
  Int memMb = round(cpu * 3.75 * 1024)
  Int disk = 100

  String output_name = output_base_name + ".fingerprintcheck"

  command <<<
    set -euo pipefail

    java -Xms${memMb - 512}m -jar /usr/gitc/picard.jar \
      CrosscheckFingerprints \
      INPUT=${write_lines(gvcf_paths)} \
      SECOND_INPUT=${write_lines(vcf_paths)} \
      H=${haplotype_database} \
      CROSSCHECK_BY=SAMPLE \
      CROSSCHECK_MODE=CHECK_SAME_SAMPLE \
      NUM_THREADS=${cpu} \
      OUTPUT=${output_name}
  >>>

  runtime {
    memory: memMb + " MB"
    disks: "local-disk " + disk + " HDD"
    preemptible: 5
  }

  output {
    File output_check = output_name
  }
}
