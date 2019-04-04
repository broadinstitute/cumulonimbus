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
    VariantsSummary phenotypeVariantsSummary
    Float pValueLimit

  }

}