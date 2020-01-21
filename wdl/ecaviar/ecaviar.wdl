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

  String regions_file_name = "regions"

  call get_regions_around_significance {
    input:
      phenotype_summary = phenotype_variants_summary.file,
      p_value_col = phenotype_variants_summary.p_value_col,
      p_value_threshold = p_value_limit,
      chromosome_col = phenotype_variants_summary.chromosome_col,
      position_col = phenotype_variants_summary.position_col,
      radius = region_padding,
      regions_file_name = regions_file_name
  }

  String region_start_col = "start"
  String region_end_col = "end"

  scatter(region in read_objects(get_regions_around_significance.regions_file)) {
    String chromosome = region[phenotype_variants_summary.chromosome_col]
    Int start = region[region_start_col]
    Int end = region[region_end_col]
    String region_notation = chromosome + "_" + start + "-" + end
    call data_munging_per_region {
      input:
        phenotype_samples_files = phenotype_samples,
        expression_samples_files = expression_data.samples_files,
        phenotype_summary = phenotype_variants_summary.file,
        id_col = phenotype_variants_summary.variant_id_col,
        chromosome_col = phenotype_variants_summary.chromosome_col,
        position_col = phenotype_variants_summary.position_col,
        phenotype_samples_name = "samples_phenotype_" + region_notation + ".vcf",
        expression_samples_name = "samples_expression_" + region_notation + ".vcf",
        chromosome = chromosome,
        start = start,
        end = end,
        phenotype_summary_sorted_name = "phenotype_summary_sorted_" + region_notation,
        common_variants_name = "common_variants" + region_notation
    }
    scatter(tissue in expression_data.tissues) {
      call clip_eqtl_region_and_get_genes {
        input:
          summary_file = tissue.summary.file,
          id_col = tissue.summary.variant_id_col,
          chromosome = chromosome,
          start = start,
          end = end,
          eqtl_region_file_name = "summary_" + tissue.tissue_name + "_" + region_notation,
          gene_id_col = tissue.gene_id_col,
          genes_file_name = "genes_" + tissue.tissue_name + "_" + region_notation
      }
      if(clip_eqtl_region_and_get_genes.n_genes > 0) {
        scatter(gene_entry in read_objects(clip_eqtl_region_and_get_genes.genes_file)) {
          String gene_id = gene_entry[tissue.gene_id_col]
          String cohort_name = gene_id + "_" + tissue.tissue_name + "_" + region_notation
          call ecaviar {
            input:
              unsorted_all2_tsv = clip_eqtl_region_and_get_genes.eqtl_region_file,
              value_col2 = tissue.gene_id_col,
              value2 = gene_id,
              position_col2 = tissue.summary.position_col,
              intersection_all_but_tsv2 = data_munging_per_region.common_variants,
              intersection_id_col = phenotype_variants_summary.variant_id_col,
              region1_tsv = data_munging_per_region.phenotype_summary_sorted,
              id_col1 = phenotype_variants_summary.variant_id_col,
              p_col1 = phenotype_variants_summary.p_value_col,
              id_col2 = tissue.summary.variant_id_col,
              p_col2 = tissue.summary.p_value_col,
              region_vcf1 = data_munging_per_region.phenotype_samples,
              region_vcf2 = data_munging_per_region.expression_samples,
              out_files_base_name = "ecaviar_" + cohort_name
          }
        }
      }
    }
  }
}

task get_regions_around_significance {
  input {
    File phenotype_summary
    String p_value_col
    Float p_value_threshold
    String chromosome_col
    String position_col
    Int radius
    String regions_file_name
  }
  runtime {
    preemptible: 3
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:200121"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    set -e
    echo "= = = $(date) = = = Now finding significant variants = = ="
    chowser tsv range --in ~{phenotype_summary} --out significant_variants.tsv \
      --col ~{p_value_col} --lt ~{p_value_threshold}
    echo "= = = $(date) = = = Beginning of significant variants = = ="
    head significant_variants.tsv
    echo "= = = $(date) = = = Now finding regions around significant variants = = ="
    chowser variants regions --in significant_variants.tsv --out ~{regions_file_name} --chrom-col ~{chromosome_col} \
     --pos-col ~{position_col} --radius ~{radius}
    echo "= = = $(date) = = = Beginning of regions = = ="
    head ~{regions_file_name}
    echo "= = = $(date) = = = Done with this task = = ="
  >>>
  output {
    File regions_file = regions_file_name
  }
}

task data_munging_per_region {
  input {
    SamplesFiles phenotype_samples_files
    SamplesFiles expression_samples_files
    File phenotype_summary
    String id_col
    String chromosome_col
    String position_col
    String chromosome
    Int start
    Int end
    String phenotype_samples_name
    String expression_samples_name
    String phenotype_summary_sorted_name
    String common_variants_name
  }
  runtime {
    preemptible: 3
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:200121"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 30 HDD"
  }
  command <<<
    set -e
    echo "= = = $(date) = = = Now clipping phenotype summary data to region = = ="
    chowser variants for-region --in ~{phenotype_summary} --out phenotype_summary_region.tsv \
      --chrom-col ~{chromosome_col} --pos-col ~{position_col} \
      --chrom ~{chromosome} --start ~{start} --end ~{end}
    echo "= = = $(date) = = = Beginning of clipped summary data = = ="
    head phenotype_summary_region.tsv
    echo "= = = $(date) = = = Now sorting phenotype summary data = = ="
    chowser tsv sort --in phenotype_summary_region.tsv --out ~{phenotype_summary_sorted_name} --col ~{position_col}
    echo "= = = $(date) = = = Sorted phenotype summary data = = ="
    head ~{phenotype_summary_sorted_name}
    echo "= = = $(date) = = = Now clipping phenotype sample data to region = = ="
    tabix -h ~{phenotype_samples_files.vcf_bgz} ~{chromosome}:~{start}-~{end-1} >~{phenotype_samples_name}
    echo "= = = $(date) = = = Now clipping expression sample data to region = = ="
    tabix -h ~{expression_samples_files.vcf_bgz} ~{chromosome}:~{start}-~{end-1} >~{expression_samples_name}
    echo "= = = $(date) = = = Now intersecting phenotype summary with phenotype samples = = ="
    chowser variants match-vcf-tsv \
      --vcf ~{phenotype_samples_name} --tsv ~{phenotype_summary_sorted_name} --id-col ~{id_col}  \
      --in-both phenotype_common_variants.tsv --vcf-only only_in_phenotype_samples \
      --tsv-only only_in_phenotype_summary
    echo "= = = $(date) = = = Beginning of intersection of phenotype summary and phenotype samples = = ="
    head phenotype_common_variants.tsv
    echo "= = = $(date) = = = Now intersecting phenotype common variants with expression samples = = ="
    chowser variants match-vcf-tsv \
      --vcf ~{expression_samples_name} --tsv phenotype_common_variants.tsv --id-col ~{id_col}  \
      --in-both ~{common_variants_name} --vcf-only only_in_expression_samples \
      --tsv-only only_in_phenotype_data
    echo "= = = $(date) = = = Beginning of intersection of phenotype variants with expression samples = = ="
    head ~{common_variants_name}
    echo "= = = $(date) = = = Done with this task = = ="
  >>>
  output {
    File phenotype_samples = phenotype_samples_name
    File expression_samples =  expression_samples_name
    File phenotype_summary_sorted = phenotype_summary_sorted_name
    File common_variants = common_variants_name
  }
}

task clip_eqtl_region_and_get_genes {
  input {
    File summary_file
    String id_col
    String chromosome
    Int start
    Int end
    String eqtl_region_file_name
    String gene_id_col
    String genes_file_name
  }
  String n_genes_file_name = "n_genes"
  runtime {
    preemptible: 3
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:200121"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    set -e
    echo "= = = $(date) = = = Clipping region from EQTL file = = ="
    chowser variants for-region-by-id --in ~{summary_file} --out ~{eqtl_region_file_name} \
      --id-col ~{id_col} --chrom ~{chromosome} --start ~{start} --end ~{end}
    echo "= = = $(date) = = = Beginning of clipped region = = ="
    head ~{eqtl_region_file_name}
    echo "= = = $(date) = = = Extracting list of genes = = ="
    chowser tsv extract-unique --in ~{eqtl_region_file_name} --out ~{genes_file_name} --col ~{gene_id_col}
    echo "= = = $(date) = = = Beginning of genes = = ="
    head ~{genes_file_name}
    echo "= = = $(date) = = = Counting number of genes = = ="
    cat ~{genes_file_name} | tail +2 | wc -l > ~{n_genes_file_name}
    echo "Number of genes is:"
    cat ~{n_genes_file_name}
    echo "= = = $(date) = = = Done with this task = = ="
  >>>
  output {
    File eqtl_region_file = eqtl_region_file_name
    File genes_file = genes_file_name
    Int n_genes = read_int(n_genes_file_name)
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
    preemptible: 3
    docker: "gcr.io/v2f-public-resources/cumulonimbus-ecaviar:200121"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    set -e
    echo "= = = $(date) = = = Filter second tsv by given value = = ="
    chowser tsv slice --in ~{unsorted_all2_tsv} --out unsorted2.tsv --col ~{value_col2} --value ~{value2}
    echo "= = = $(date) = = = Beginning of filtered second TSV = = ="
    head unsorted2.tsv
    echo "= = = $(date) = = = Sort second tsv by position = = ="
    chowser tsv sort-ids --in unsorted2.tsv --out region2.tsv --col ~{id_col2}
    echo "= = = $(date) = = = Beginning of sorted second TSV = = ="
    head region2.tsv
    echo "= = = $(date) = = = Beginning of intersection of all but second TSV = = ="
    head ~{intersection_all_but_tsv2}
    echo "= = = $(date) = = = Intersect common variants with second tsv = = ="
    chowser variants match-tsv-tsv \
      --tsv1 ~{intersection_all_but_tsv2} --tsv2 region2.tsv --id-col1 ~{intersection_id_col} --id-col2 ~{id_col2}  \
      --in-both selected.tsv --tsv1-only unmatched1.tsv --tsv2-only unmatched2.tsv
    echo "= = = $(date) = = = Beginning of selected variants = = ="
    head selected.tsv
    echo "= = = $(date) = = = Checking if there is more than one variant = = ="
    if [ $(tail +2 selected.tsv|wc -l) -gt 1 ]; then
      echo "= = = $(date) = = = There is more than one variant = = ="
      echo "= = = $(date) = = = Filter selected variants from first region (tsv) = = ="
      chowser variants select-tsv --data ~{region1_tsv} --selection selected.tsv --out selected1.tsv \
      --id-col-data ~{id_col1} --id-col-selection ~{intersection_id_col}
      echo "= = = $(date) = = = Beginning of selected variants from first region (tsv) = = ="
      head selected1.tsv
      echo "= = = $(date) = = = Filter selected variants from second region (tsv) = = ="
      chowser variants select-tsv --data region2.tsv --selection selected.tsv --out selected2.tsv \
      --id-col-data ~{id_col2} --id-col-selection ~{intersection_id_col}
      echo "= = = $(date) = = = Beginning of selected variants from second region (tsv) = = ="
      head selected2.tsv
      echo "= = = $(date) = = = Filter selected variants from first region (vcf) = = ="
      chowser variants select-vcf --data ~{region_vcf1} --selection selected.tsv --out selected1.vcf \
        --id-col-selection ~{intersection_id_col}
      echo "= = = $(date) = = = Beginning of selected variants from first region (vcf) = = ="
      head -n 20 selected1.vcf
      echo "= = = $(date) = = = Filter selected variants from second region (vcf) = = ="
      chowser variants select-vcf --data ~{region_vcf2} --selection selected.tsv --out selected2.vcf \
        --id-col-selection ~{intersection_id_col}
      echo "= = = $(date) = = = Beginning of selected variants from second region (vcf) = = ="
      head -n 20 selected2.vcf
      echo "= = = $(date) = = = Calculating first set of correlations = = ="
      plink --vcf selected1.vcf --a1-allele selected1.vcf 4 3 '#' --r --out plink1
      echo "= = = $(date) = = = First set of correlations = = ="
      head plink1.ld
      echo "= = = $(date) = = = Casting first set of correlations into a square matrix = = ="
      chowser caviar matrix --ids-file selected1.vcf --values-file plink1.ld \
        --value-col R --id-col1 SNP_A --id-col2 SNP_B --out ld_file1
      echo "= = = $(date) = = = First set of correlations as square matrix= = ="
      head ld_file1
      echo "= = = $(date) = = = Calculating second set of correlations = = ="
      plink --vcf selected2.vcf --a1-allele selected2.vcf 4 3 '#' --r --out plink2
      echo "= = = $(date) = = = Second set of correlations = = ="
      head plink2.ld
      echo "= = = $(date) = = = Casting second set of correlations into a square matrix = = ="
      chowser caviar matrix --ids-file selected2.vcf --values-file plink2.ld \
        --value-col R --id-col1 SNP_A --id-col2 SNP_B --out ld_file2
      echo "= = = $(date) = = = Second set of correlations as square matrix= = ="
      head ld_file2
      echo "= = = $(date) = = = Calculating first file of Z-values = = ="
      chowser caviar p-to-z --in selected1.tsv --out z_file1 --id-col ~{id_col1} --p-col ~{p_col1}
      echo "= = = $(date) = = = First file of Z-values = = ="
      head z_file1
      echo "= = = $(date) = = = Calculating second file of Z-values = = ="
      chowser caviar p-to-z --in selected2.tsv --out z_file2 --id-col ~{id_col2} --p-col ~{p_col2}
      echo "= = = $(date) = = = Second file of Z-values = = ="
      head z_file2
      echo "= = = $(date) = = = Running eCAVIAR = = ="
      eCAVIAR -l ld_file1 -z z_file1 -l ld_file2 -z z_file2 -o ~{out_files_base_name}
    else
      echo "= = = $(date) = = = There is not more than one variant  = = ="
    fi
    echo "= = = $(date) = = = Done with this task = = ="
  >>>
  output {
    Array[File] out_file_col = glob(out_files_base_name + "_col")
    Array[File] out_file_1_set = glob(out_files_base_name + "_1_set")
    Array[File] out_file_1_post = glob(out_files_base_name + "_1_post")
    Array[File] out_file_2_set = glob(out_files_base_name + "_2_set")
    Array[File] out_file_2_post = glob(out_files_base_name + "_2_post")
  }
}
