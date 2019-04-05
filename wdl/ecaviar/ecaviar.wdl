version 1.0

struct SamplesFiles {
  File VcfGz
  File Tbi
}

struct VariantsSummary {
  File file
  String variant_id_col
  String p_value_col
  String chromosome_col
  String position_col
}

workflow ecaviar {
  input {
    SamplesFiles phenotype_samples
    VariantsSummary phenotype_variants_summary
    Float p_value_limit
    Int region_padding

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
    call region_test {
      input:
        chromosome = chromosome,
        start = start,
        end = end,
        out_file_name = "region"
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
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190403"
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
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190403"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    chowser variants regions --in ~{in_file} --out ~{out_file_name} --chrom ~{chromosome_col} \
     --pos ~{position_col} --radius ~{radius}
  >>>
  output {
    File out_file = out_file_name
  }
}

task region_test {
  input {
    String chromosome
    Int start
    Int end
    String out_file_name
  }
  runtime {
    docker: "gcr.io/broad-gdr-dig-storage/cumulonimbus-ecaviar:190403"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    echo "Chromosome: ~{chromosome}, start: ~{start}, end: ~{end}" > ~{out_file_name}
  >>>
  output {
    File out_file = out_file_name
  }
}