version 1.0

struct SamplesFiles {
  File vcf_bgz
  File vcf_bgz_tbi
}

struct VariantsSummary {
  File file
  String variant_id_col
  String p_value_col
  String chromosome_col
  String position_col
}

struct Tissue {
  File variants
}

workflow ecaviar {
  input {
    SamplesFiles phenotype_samples
    VariantsSummary phenotype_variants_summary
    Float p_value_limit
    Int region_padding
    Array[Tissue] tissues
  }

  String significant_variants_file_name = "significant_variants"

  call get_phenotype_significant_variants {
    input:
      in_file = phenotype_variants_summary.file,
      p_value_col = phenotype_variants_summary.p_value_col,
      p_value_limit = p_value_limit,
      out_file_name = significant_variants_file_name
  }

  String regions_file_name = "regions"

  call get_regions_around_significance {
    input:
      in_file = get_phenotype_significant_variants.out_file,
      chromosome_col = phenotype_variants_summary.chromosome_col,
      position_col = phenotype_variants_summary.position_col,
      radius = region_padding,
      out_file_name = regions_file_name
  }

  String regionStartCol = "start"
  String regionEndCol = "end"

  scatter(region in read_objects(get_regions_around_significance.out_file)) {
    String chromosome = region[phenotype_variants_summary.chromosome_col]
    Int start = region[regionStartCol]
    Int end = region[regionEndCol]
    call clip_region_from_samples as clip_region_from_phenotype_samples {
      input:
        samples_files = phenotype_samples,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "samples_" + chromosome + ":" + start + "-" + end + ".vcf"
    }
    call clip_region_from_summary as clip_region_from_phenotype_summary {
      input:
        in_file = get_phenotype_significant_variants.out_file,
        chromosome_col = phenotype_variants_summary.chromosome_col,
        position_col = phenotype_variants_summary.position_col,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "summary_" + chromosome + ":" + start + "-" + end + ".tsv"
    }
    call sort_file_by_col {
      input:
        in_file = clip_region_from_phenotype_summary.out_file,
        col = phenotype_variants_summary.position_col,
        out_file_name = "summary_sorted_" + chromosome + ":" + start + "-" + end + ".tsv"
    }
    call match_variants {
      input:
        in_vcf = clip_region_from_phenotype_samples.out_file,
        in_tsv = clip_region_from_phenotype_summary.out_file,
        id_col = phenotype_variants_summary.variant_id_col,
        out_both_name = "variants_common_" + chromosome + ":" + start + "-" + end + ".tsv",
        out_vcf_only_name = "variants_samples_only_" + chromosome + ":" + start + "-" + end + ".tsv",
        out_tsv_only_name = "variants_summary_only_" + chromosome + ":" + start + "-" + end + ".tsv"
    }
    scatter(tissue in tissues) {
      call clip_region_from_summary as clip_region_from_expression_summary {
        input:
          in_file = tissue.variants
      }
    }
  }
}

task get_phenotype_significant_variants {
  input {
    File in_file
    String p_value_col
    Float p_value_limit
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190424"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser tsv filter --in ~{in_file} --out ~{out_file_name} --col ~{p_value_col} --lt ~{p_value_limit}
  >>>
  output {
    File out_file = out_file_name
  }
}

task get_regions_around_significance {
  input {
    File in_file
    String chromosome_col
    String position_col
    Int radius
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190424"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants regions --in ~{in_file} --out ~{out_file_name} --chrom-col ~{chromosome_col} \
     --pos-col ~{position_col} --radius ~{radius}
  >>>
  output {
    File out_file = out_file_name
  }
}

task clip_region_from_samples {
  input {
    SamplesFiles samples_files
    String chromosome
    Int start
    Int end
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190424"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    tabix -h ~{samples_files.vcf_bgz} ~{chromosome}:~{start}-~{end-1} >~{out_file_name}
  >>>
  output {
    File out_file = out_file_name
  }
}

task clip_region_from_summary {
  input {
    File in_file
    String chromosome_col
    String position_col
    String out_file_name
    String chromosome
    Int start
    Int end
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190424"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants for-region --in ~{in_file} --out ~{out_file_name} \
      --chrom-col ~{chromosome_col} --pos-col ~{position_col} --chrom ~{chromosome} --start ~{start} --end ~{end}
  >>>
  output {
    File out_file = out_file_name
  }
}

task sort_file_by_col {
  input {
    File in_file
    String col
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190424"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser tsv sort --in ~{in_file} --out ~{out_file_name} --col ~{col}
  >>>
  output {
    File out_file = out_file_name
  }
}

task match_variants {
  input {
    File in_vcf
    File in_tsv
    String id_col
    String out_both_name
    String out_vcf_only_name
    String out_tsv_only_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190424"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser match variants \
      --vcf ~{in_vcf} --tsv ~{in_tsv} --id-col ~{id_col}  \
      --in-both ~{out_both_name} --vcf-only ~{out_vcf_only_name} --tsv-only ~{out_tsv_only_name}
  >>>
  output {
    File out_both = out_both_name
    File out_vcf_only = out_vcf_only_name
    File out_tsv_only = out_tsv_only_name
  }
}