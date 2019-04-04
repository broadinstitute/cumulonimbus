version 1.0

struct SamplesFiles {
  File VcfGz
  File Tbi
}

struct VariantsSummary {
  File file
  String variantIdCol
  String pValueCol
  String chromosomeCol
  String positionCol
}

workflow ecaviar {
  input {
    SamplesFiles phenotypeSamples
    VariantsSummary phenotype_variants_summary
    Float pValueLimit

  }

  call get_phenotype_significant_variants {
    input:
      variants_summary = phenotype_variants_summary
  }
}

task get_phenotype_significant_variants {
  input {
    VariantsSummary variants_summary
    String out_file_name
  }
  command <<<
    chowser tsv filter
  >>>
}