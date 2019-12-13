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
  String ref_col
  String alt_col
}

struct Tissue {
  String tissue_name
  String gene_id_col
  VariantsSummary summary
}

struct ExpressionData {
  SamplesFiles samples_files
  Array[Tissue] tissues
}

workflow ecaviar {
  input {
    SamplesFiles phenotype_samples
    VariantsSummary phenotype_variants_summary
    Float p_value_limit
    Int region_padding
    ExpressionData expression_data
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

  String region_start_col = "start"
  String region_end_col = "end"

  scatter(region in read_objects(get_regions_around_significance.out_file)) {
    String chromosome = region[phenotype_variants_summary.chromosome_col]
    Int start = region[region_start_col]
    Int end = region[region_end_col]
    String region_notation = chromosome + "_" + start + "-" + end
    call clip_region_from_samples as clip_region_from_phenotype_samples {
      input:
        samples_files = phenotype_samples,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "samples_phenotype_" + region_notation + ".vcf"
    }
    call clip_region_from_samples as clip_region_from_expression_samples {
      input:
        samples_files = expression_data.samples_files,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "samples_expression_" + region_notation + ".vcf"
    }
    call clip_region_from_summary as clip_region_from_phenotype_summary {
      input:
        in_file = get_phenotype_significant_variants.out_file,
        chromosome_col = phenotype_variants_summary.chromosome_col,
        position_col = phenotype_variants_summary.position_col,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "summary_" + region_notation
    }
    call sort_file_by_col as sort_phenotype_summary {
      input:
        in_file = clip_region_from_phenotype_summary.out_file,
        col = phenotype_variants_summary.position_col,
        out_file_name = "phenotype_summary_sorted_" + region_notation
    }
    call match_variants_vcf_tsv as match_variants_phenotype_samples_summary {
      input:
        in_vcf = clip_region_from_phenotype_samples.out_file,
        in_tsv = sort_phenotype_summary.out_file,
        id_col = phenotype_variants_summary.variant_id_col,
        out_both_name = "phenotype_common_" + region_notation,
        out_vcf_only_name = "phenotype_samples_only_" + region_notation,
        out_tsv_only_name = "phenotype_summary_only_" + region_notation
    }
    call match_variants_vcf_tsv as match_variants_phenotype_expression_samples {
      input:
        in_vcf = clip_region_from_expression_samples.out_file,
        in_tsv = match_variants_phenotype_samples_summary.out_both,
        id_col = phenotype_variants_summary.variant_id_col,
        out_both_name = "phenotype_expression_samples_" + region_notation,
        out_vcf_only_name = "expression_samples_only_" + region_notation,
        out_tsv_only_name = "phenotype_only_" + region_notation
    }
    scatter(tissue in expression_data.tissues) {
      call clip_region_from_summary as clip_region_from_expression_summary {
        input:
          in_file = tissue.summary.file,
          chromosome_col = tissue.summary.chromosome_col,
          position_col = tissue.summary.position_col,
          chromosome = chromosome,
          start = start,
          end = end,
          out_file_name = "summary_" + tissue.tissue_name + "_" + region_notation
      }
      Int n_expression_entries = clip_region_from_expression_summary.count
      String report_by_region_tissue =
        "Tissue " + tissue.tissue_name + " in region " + region_notation + " has " +
        n_expression_entries + " expression entries."
      if(n_expression_entries > 0) {
        call extract_unique as extract_genes {
          input:
            in_file = clip_region_from_expression_summary.out_file,
            col = tissue.gene_id_col,
            out_file_name = "genes_" + tissue.tissue_name + "_" + region_notation
        }
        scatter(gene_entry in read_objects(extract_genes.out_file)) {
          String gene_id = gene_entry[tissue.gene_id_col]
          String cohort_name = gene_id + "_" + tissue.tissue_name + "_" + region_notation
          call ecaviar {
            input:
              unsorted_all2_tsv = clip_region_from_expression_summary.out_file,
              value_col2 = tissue.gene_id_col,
              value2 = gene_id,
              position_col2 = tissue.summary.position_col,
              intersection_all_but_tsv2 = match_variants_phenotype_expression_samples.out_both,
              intersection_id_col = tissue.summary.variant_id_col,
              region1_tsv = sort_phenotype_summary.out_file,
              id_col1 = phenotype_variants_summary.variant_id_col,
              p_col1 = phenotype_variants_summary.p_value_col,
              id_col2 = tissue.summary.variant_id_col,
              p_col2 = tissue.summary.p_value_col,
              region_vcf1 = clip_region_from_phenotype_samples.out_file,
              region_vcf2 = clip_region_from_expression_samples.out_file,
              out_files_base_name = "ecaviar_" + cohort_name
          }
        }
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
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser tsv range --in ~{in_file} --out ~{out_file_name} --col ~{p_value_col} --lt ~{p_value_limit}
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
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
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
  String count_suffix = ".count"
  String out_count_file_name = out_file_name + count_suffix
  runtime {
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    tabix -h ~{samples_files.vcf_bgz} ~{chromosome}:~{start}-~{end-1} >~{out_file_name}
    grep "#" -v ~{out_file_name} |wc -l > ~{out_count_file_name}
  >>>
  output {
    File out_file = out_file_name
    Int count = read_int(out_count_file_name)
  }
}

task clip_region_from_summary {
  input {
    File in_file
    String out_file_name
    String chromosome_col
    String position_col
    String chromosome
    Int start
    Int end
  }
  String count_suffix = ".count"
  String out_count_file_name = out_file_name + count_suffix
  runtime {
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants for-region --in ~{in_file} --out ~{out_file_name} \
      --chrom-col ~{chromosome_col} --pos-col ~{position_col} \
      --chrom ~{chromosome} --start ~{start} --end ~{end}
    tail -n +2 ~{out_file_name} | wc -l > ~{out_count_file_name}
  >>>
  output {
    File out_file = out_file_name
    Int count = read_int(out_count_file_name)
  }
}

task sort_file_by_col {
  input {
    File in_file
    String col
    String out_file_name
  }
  runtime {
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
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

task match_variants_vcf_tsv {
  input {
    File in_vcf
    File in_tsv
    String id_col
    String out_both_name
    String out_vcf_only_name
    String out_tsv_only_name
  }
  String count_suffix = ".count"
  String out_both_count_name = out_both_name + count_suffix
  String out_vcf_only_count_name = out_vcf_only_name + count_suffix
  String out_tsv_only_count_name = out_tsv_only_name + count_suffix
  runtime {
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants match-vcf-tsv \
      --vcf ~{in_vcf} --tsv ~{in_tsv} --id-col ~{id_col}  \
      --in-both ~{out_both_name} --vcf-only ~{out_vcf_only_name} --tsv-only ~{out_tsv_only_name}
    tail -n +2 ~{out_both_name} | wc -l > ~{out_both_count_name}
    tail -n +2 ~{out_vcf_only_name} | wc -l > ~{out_vcf_only_count_name}
    tail -n +2 ~{out_tsv_only_name} | wc -l > ~{out_tsv_only_count_name}
  >>>
  output {
    File out_both = out_both_name
    Int out_both_count = read_int(out_both_count_name)
    File out_vcf_only = out_vcf_only_name
    Int out_vcf_only_count = read_int(out_vcf_only_count_name)
    File out_tsv_only = out_tsv_only_name
    Int out_tsv_only_count = read_int(out_tsv_only_count_name)
  }
}

task extract_unique {
  input {
    File in_file
    String col
    String out_file_name
  }
  runtime {
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser tsv extract-unique --in ~{in_file} --out ~{out_file_name} --col ~{col}
  >>>
  output {
    File out_file = out_file_name
  }
}

task ecaviar {
  input {
    File unsorted_all2_tsv
    String value_col2
    String value2
    String position_col2
    File intersection_all_but_tsv2
    String intersection_id_col
    File region1_tsv
    String id_col1
    String p_col1
    String id_col2
    String p_col2
    File region_vcf1
    File region_vcf2
    String out_files_base_name
  }
  runtime {
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:191206"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    echo "= = = Filter second all tsv by given value = = ="
    chowser tsv slice --in ~{unsorted_all2_tsv} --out unsorted2.tsv --col ~{value_col2} --value ~{value2}
    echo "= = = Sort second tsv by position = = ="
    chowser tsv sort --in unsorted2.tsv --out region2.tsv --col ~{position_col2}
    echo "= = = Intersect common variants with second tsv = = ="
    chowser variants match-tsv-tsv \
      --tsv1 region2.tsv --tsv2 ~{intersection_all_but_tsv2} --id-col1 ~{id_col2} --id-col2 ~{intersection_id_col}  \
      --in-both selected.tsv --tsv1-only unmatched1.tsv --tsv2-only unmatched2.tsv
    echo "= = = Checking if there is more than one variant = = ="
    if [ $(tail +2 selected.tsv|wc -l) -gt 1 ]; then
      echo "= = = There is more than one variant = = ="
      echo "= = = Filter selected variants from first region (tsv) = = ="
      chowser variants select-tsv --data ~{region1_tsv} --selection selected.tsv --out selected1.tsv \
      --id-col-data ~{id_col1} --id-col-selection ~{intersection_id_col}
      echo "= = = Filter selected variants from second region (tsv) = = ="
      chowser variants select-tsv --data region2.tsv --selection selected.tsv --out selected2.tsv \
      --id-col-data ~{id_col2} --id-col-selection ~{intersection_id_col}
      echo "= = = Filter selected variants from first region (vcf) = = ="
      chowser variants select-vcf --data ~{region_vcf1} --selection selected.tsv --out selected1.vcf \
        --id-col-selection ~{intersection_id_col}
      echo "= = = Filter selected variants from second region (vcf) = = ="
      chowser variants select-vcf --data ~{region_vcf2} --selection selected.tsv --out selected2.vcf \
        --id-col-selection ~{intersection_id_col}
      echo "= = = Calculating first set of correlations = = ="
      plink --vcf selected1.vcf --a1-allele selected1.vcf 4 3 '#' --r --out plink
      echo "= = = Casting first set of correlations into a square matrix = = ="
      chowser caviar matrix --ids-file selected1.vcf --values-file plink.ld \
        --value-col R --id-col1 SNP_A --id-col2 SNP_B --out ld_file1
      echo "= = = Calculating second set of correlations = = ="
      plink --vcf selected2.vcf --a1-allele selected2.vcf 4 3 '#' --r --out plink
      echo "= = = Casting second set of correlations into a square matrix = = ="
      chowser caviar matrix --ids-file selected2.vcf --values-file plink.ld \
        --value-col R --id-col1 SNP_A --id-col2 SNP_B --out ld_file2
      echo "= = = Calculating first file of Z-values = = ="
      chowser caviar p-to-z --in selected1.tsv --out z_file1 --id-col ~{id_col1} --p-col ~{p_col1}
      echo "= = = Calculating second file of Z-values = = ="
      chowser caviar p-to-z --in selected2.tsv --out z_file2 --id-col ~{id_col2} --p-col ~{p_col2}
      echo "= = = Running eCAVIAR = = ="
      eCAVIAR -l ld_file1 -z z_file1 -l ld_file2 -z z_file2 -o ~{out_files_base_name}
    else
      echo "= = = There is not more than one variant  = = ="
    fi
    echo "= = = Done with this task = = ="
  >>>
  output {
    File? out_file_col = out_files_base_name + "_col"
    File? out_file_1_set = out_files_base_name + "_1_set"
    File? out_file_1_post = out_files_base_name + "_1_post"
    File? out_file_2_set = out_files_base_name + "_2_set"
    File? out_file_2_post = out_files_base_name + "_2_post"
  }
}
