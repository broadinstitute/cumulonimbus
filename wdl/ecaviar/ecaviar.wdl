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

  String canonicalized_phenotype_samples_file_base_name = "samples_phenotype_canon"

  call canonicalize_samples as canonicalize_phenotype_samples {
    input:
      in_files = phenotype_samples,
      out_files_base_name = canonicalized_phenotype_samples_file_base_name
  }

  String canonicalized_expression_samples_file_base_name = "samples_expression_canon"

  call canonicalize_samples as canonicalize_expression_samples {
    input:
      in_files = expression_data.samples_files,
      out_files_base_name = canonicalized_expression_samples_file_base_name
  }

  String phenotype_summary_canonicalized_name = "phenotype_summary_canon"

  call canonicalize_summary as canonicalize_phenotype_summary {
    input:
      summary = phenotype_variants_summary,
      out_file_name = phenotype_summary_canonicalized_name
  }

  String significant_variants_file_name = "significant_variants"

  call get_phenotype_significant_variants {
    input:
      in_file = canonicalize_phenotype_summary.summary_canonical.file,
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
        samples_files = canonicalize_phenotype_samples.out_files,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "samples_phenotype_" + chromosome + ":" + start + "-" + end + ".vcf"
    }
    call clip_region_from_samples as clip_region_from_expression_samples {
      input:
        samples_files = canonicalize_expression_samples.out_files,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "samples_expression_" + chromosome + ":" + start + "-" + end + ".vcf"
    }
    call clip_region_from_summary as clip_region_from_phenotype_summary {
      input:
        in_file = get_phenotype_significant_variants.out_file,
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "summary_" + chromosome + ":" + start + "-" + end
    }
    call sort_file_by_col as sort_phenotype_summary {
      input:
        in_file = clip_region_from_phenotype_summary.out_file,
        col = phenotype_variants_summary.position_col,
        out_file_name = "phenotype_summary_sorted_" + chromosome + ":" + start + "-" + end
    }
    call match_variants_vcf_tsv as match_variants_phenotype_samples_summary {
      input:
        in_vcf = clip_region_from_phenotype_samples.out_file,
        in_tsv = sort_phenotype_summary.out_file,
        id_col = phenotype_variants_summary.variant_id_col,
        out_both_name = "phenotype_common_" + chromosome + ":" + start + "-" + end,
        out_vcf_only_name = "phenotype_samples_only_" + chromosome + ":" + start + "-" + end,
        out_tsv_only_name = "phenotype_summary_only_" + chromosome + ":" + start + "-" + end
    }
    call match_variants_vcf_tsv as match_variants_phenotype_expression_samples {
      input:
        in_vcf = clip_region_from_expression_samples.out_file,
        in_tsv = match_variants_phenotype_samples_summary.out_both,
        id_col = phenotype_variants_summary.variant_id_col,
        out_both_name = "phenotype_expression_samples_" + chromosome + ":" + start + "-" + end,
        out_vcf_only_name = "expression_samples_only_" + chromosome + ":" + start + "-" + end,
        out_tsv_only_name = "phenotype_only_" + chromosome + ":" + start + "-" + end
    }
    scatter(tissue in expression_data.tissues) {
      call clip_region_from_summary as clip_region_from_expression_summary {
        input:
          in_file = tissue.summary.file,
          chromosome = chromosome,
          start = start,
          end = end,
          out_file_name = "summary_" + tissue.tissue_name + "_" + chromosome + ":" + start + "-" + end
      }
      call extract_unique as extract_genes {
        input:
          in_file = clip_region_from_expression_summary.out_file,
          col = tissue.gene_id_col,
          out_file_name = "genes_" + tissue.tissue_name + "_" + chromosome + ":" + start + "-" + end
      }
      scatter(gene_entry in read_objects(extract_genes.out_file)) {
        String gene_id = gene_entry[tissue.gene_id_col]
        String cohort_name = gene_id + "_" + tissue.tissue_name + "_" + chromosome + ":" + start + "-" + end
        call slice_by_value as slice_by_gene_id {
          input:
            in_file = clip_region_from_expression_summary.out_file,
            col = tissue.gene_id_col,
            value = gene_id,
            out_file_name = "summary_" + cohort_name
        }
        call sort_file_by_col as sort_cohort_by_position {
          input:
            in_file = slice_by_gene_id.out_file,
            col = tissue.summary.position_col,
            out_file_name = "summary_sorted_" + cohort_name
        }
        call match_variants_tsv_tsv as match_variants_with_expression_summary {
          input:
            in_tsv1 = match_variants_phenotype_expression_samples.out_both,
            in_tsv2 = sort_cohort_by_position.out_file,
            id_col1 = tissue.summary.variant_id_col,
            id_col2 = tissue.summary.variant_id_col,
            out_both_name = "variants_" + cohort_name,
            out_tsv1_only_name = "variants_not_in_expression_summary_" + cohort_name,
            out_tsv2_only_name = "variants_only_in_expression_summary_" + cohort_name,
        }
        call select_variants_tsv as select_variants_phenotype_summary {
          input:
            data_file = sort_phenotype_summary.out_file,
            selection_file = match_variants_with_expression_summary.out_both,
            col_in_data = phenotype_variants_summary.variant_id_col,
            col_in_selection = tissue.summary.variant_id_col,
            out_file_name = "phenotype_summary_selected_variants_" + cohort_name
        }
        call select_variants_tsv as select_variants_expression_summary {
          input:
            data_file = sort_cohort_by_position.out_file,
            selection_file = match_variants_with_expression_summary.out_both,
            col_in_data = tissue.summary.variant_id_col,
            col_in_selection = tissue.summary.variant_id_col,
            out_file_name = "expression_summary_selected_variants_" + cohort_name
        }
        call select_variants_vcf as select_variants_phenotype_samples {
          input:
            data_file = clip_region_from_phenotype_samples.out_file,
            selection_file = match_variants_with_expression_summary.out_both,
            col_in_selection = tissue.summary.variant_id_col,
            out_file_name = "phenotype_samples_selected_variants_" + cohort_name + ".vcf"
        }
        call select_variants_vcf as select_variants_expression_samples {
          input:
            data_file = clip_region_from_expression_samples.out_file,
            selection_file = match_variants_with_expression_summary.out_both,
            col_in_selection = tissue.summary.variant_id_col,
            out_file_name = "expression_samples_selected_variants_" + cohort_name + ".vcf"
        }
        call calculate_correlations as calculate_phenotype_correlations {
          input:
            in_file = select_variants_phenotype_samples.out_file,
            out_file_name = "phenotype_correlations_" + cohort_name
        }
        call calculate_correlations as calculate_expression_correlations {
          input:
            in_file = select_variants_expression_samples.out_file,
            out_file_name = "expression_correlations_" + cohort_name
        }
        call generate_z_scores_for_ecaviar as generate_z_scores_phenotype {
          input:
            in_file = select_variants_phenotype_summary.out_file,
            id_col = phenotype_variants_summary.variant_id_col,
            p_col = phenotype_variants_summary.p_value_col,
            out_file_name = "z_scores_phenotype_" + cohort_name
        }
        call generate_z_scores_for_ecaviar as generate_z_scores_expression {
          input:
            in_file = select_variants_phenotype_summary.out_file,
            id_col = tissue.summary.variant_id_col,
            p_col = tissue.summary.p_value_col,
            out_file_name = "z_scores_phenotype_" + cohort_name
        }
        call ecaviar {
          input:
            ld_file1 = calculate_phenotype_correlations.out_file,
            z_file1 = generate_z_scores_phenotype.out_file,
            ld_file2 = calculate_expression_correlations.out_file,
            z_file2 = generate_z_scores_expression.out_file,
            out_file_name = "ecaviar_" + cohort_name
        }
      }
    }
  }
}

task canonicalize_samples {
  input {
    SamplesFiles in_files
    String out_files_base_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    gunzip -c ~{in_files.vcf_bgz} >variants_temp.vcf
    chowser variants canonicalize-vcf --in variants_temp.vcf --out ~{out_files_base_name}.vcf
    bgzip -c ~{out_files_base_name}.vcf > ~{out_files_base_name}.vcf.bgz
    tabix -p vcf ~{out_files_base_name}.vcf.bgz
  >>>
  output {
    SamplesFiles out_files = {
      "vcf_bgz": out_files_base_name + ".vcf.bgz",
      "vcf_bgz_tbi": out_files_base_name + ".vcf.bgz.tbi"
    }
  }
}

task canonicalize_summary {
  input {
    VariantsSummary summary
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants canonicalize-tsv --in ~{summary.file} --out ~{out_file_name} --id-col ~{summary.variant_id_col} \
      --chrom-col ~{summary.chromosome_col} --pos-col ~{summary.position_col} --ref-col ~{summary.ref_col} \
      --alt-col ~{summary.alt_col}
  >>>
  output {
    VariantsSummary summary_canonical = {
        "file" : out_file_name,
        "variant_id_col" : summary.variant_id_col,
        "p_value_col" : summary.p_value_col,
        "chromosome_col" : summary.chromosome_col,
        "position_col" : summary.position_col,
        "ref_col" : summary.ref_col,
        "alt_col" : summary.alt_col
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
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
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
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
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
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
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
    String out_file_name
    String chromosome
    Int start
    Int end
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants for-region --in ~{in_file} --out ~{out_file_name} \
      --chrom ~{chromosome} --start ~{start} --end ~{end}
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
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
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
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants match-vcf-tsv \
      --vcf ~{in_vcf} --tsv ~{in_tsv} --id-col ~{id_col}  \
      --in-both ~{out_both_name} --vcf-only ~{out_vcf_only_name} --tsv-only ~{out_tsv_only_name}
  >>>
  output {
    File out_both = out_both_name
    File out_vcf_only = out_vcf_only_name
    File out_tsv_only = out_tsv_only_name
  }
}

task match_variants_tsv_tsv {
  input {
    File in_tsv1
    File in_tsv2
    String id_col1
    String id_col2
    String out_both_name
    String out_tsv1_only_name
    String out_tsv2_only_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants match-tsv-tsv \
      --tsv1 ~{in_tsv1} --tsv2 ~{in_tsv2} --id-col1 ~{id_col1} --id-col2 ~{id_col2}  \
      --in-both ~{out_both_name} --tsv1-only ~{out_tsv1_only_name} --tsv2-only ~{out_tsv2_only_name}
  >>>
  output {
    File out_both = out_both_name
    File out_tsv1_only = out_tsv1_only_name
    File out_tsv2_only = out_tsv2_only_name
  }
}

task extract_unique {
  input {
    File in_file
    String col
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
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

task slice_by_value {
  input {
    File in_file
    String col
    String value
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser tsv slice --in ~{in_file} --out ~{out_file_name} --col ~{col} --value ~{value}
  >>>
  output {
    File out_file = out_file_name
  }
}

task select_variants_tsv {
  input {
    File data_file
    File selection_file
    String col_in_data
    String col_in_selection
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants select-tsv --data ~{data_file} --selection ~{selection_file} --out ~{out_file_name} \
    --id-col-data ~{col_in_data} --id-col-selection ~{col_in_selection}
  >>>
  output {
    File out_file = out_file_name
  }
}

task select_variants_vcf {
  input {
    File data_file
    File selection_file
    String col_in_selection
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    bgzip ~{data_file}
    chowser variants select-vcf --data ~{data_file + ".bgz"} --selection ~{selection_file} --out ~{out_file_name} \
      --id-col-selection ~{col_in_selection}
    gunzip ~{out_file_name + ".bgz"}
  >>>
  output {
    File out_file = out_file_name
  }
}

task calculate_correlations {
  input {
    File in_file
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    plink --vcf ~{in_file} --a1-allele ~{in_file} 4 3 '#' --r --out plink
    chowser caviar matrix --ids-file ~{in_file} --values-file plink.ld \
      --value-col R --id-col1 SNP_A --id-col2 SNP_B --out ~{out_file_name}
  >>>
  output {
    File out_file = out_file_name
  }
}

task generate_z_scores_for_ecaviar {
  input {
    File in_file
    String id_col
    String p_col
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser caviar p-to-z --in ~{in_file} --out ~{out_file_name} --id-col ~{id_col} --p-col ~{p_col}
  >>>
  output {
    File out_file = out_file_name
  }
}

task ecaviar {
  input {
    File ld_file1
    File z_file1
    File ld_file2
    File z_file2
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190524"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    eCAVIAR -l ~{ld_file1} -z ~{z_file1} -l ~{ld_file2} -z ~{z_file2} -o ~{out_file_name}
  >>>
  output {
    File out_file = out_file_name
  }
}

