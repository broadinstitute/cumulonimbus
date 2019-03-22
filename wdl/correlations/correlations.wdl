version 1.0

workflow correlations {
  input {
    File variants
    File index
    String chromosome
    Int regionStart
    Int regionEnd
    String out_file_base_name
  }
  String region_file_name = "region.vcf"

  call select_region {
    input:
      variants = variants,
      index = index,
      chromosome = chromosome,
      regionStart = regionStart,
      regionEnd = regionEnd,
      out_file_name = region_file_name
  }

  call calculate_correlations {
    input:
      variants = select_region.out_file,
      out_file_base_name = out_file_base_name
  }
}

task select_region {
  input {
    File variants
    File index
    String chromosome
    Int regionStart
    Int regionEnd
    String out_file_name
  }
  runtime {
    docker: "us.gcr.io/broad-gdr-dig-storage/tabix-plink:2019-03-01"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    tabix -h ~{variants} ~{chromosome}:~{regionStart}-~{regionEnd} >~{out_file_name}
  >>>
  output {
    File out_file = out_file_name
  }
}

task calculate_correlations {
  input {
    File variants
    String out_file_base_name
  }
  runtime {
    docker: "us.gcr.io/broad-gdr-dig-storage/tabix-plink:2019-03-01"
    cpu: 1
    memory: "5 GB"
    disks: "local-disk 20 HDD"
  }
  command <<<
    plink --vcf ~{variants} --r --matrix --out ~{out_file_base_name}
  >>>
  output {
    File out_file = out_file_base_name + ".ld"
  }
}