workflow numetal {
  Array[File] inputs
  String out_prefix
  String out_postfix

  call metal {
    input:
      inputs = inputs,
      out_prefix = out_prefix,
      out_postfix = out_postfix
  }

  output {
    File out = metal.out
  }
}

task metal {
  Array[File] inputs
  String out_prefix
  String out_postfix
  runtime {
    docker: "us.gcr.io/broad-gdr-dig-storage/metal-python:2018-08-28"
    cpu: 1
    memory: "3 GB"
    disks: "local-disk 20 HDD"
  }
  command {
    cat << EOF >metalcast.py
    inputs = ["${sep='\", \"' inputs}"]
    preamble = """
    SCHEME STDERR
    MARKER ID
    PVALUE P
    FREQ MAF
    ALLELE ALT REF
    AVERAGEFREQ ON
    MINMAXFREQ ON
    STDERR SEBETA
    EFFECT BETA
    CUSTOMVARIABLE TotalSampleSize
    LABEL TotalSampleSize as NS
    """
    processes = reduce(lambda item1, item2: item1 + item2, map(lambda item: "PROCESS " + item + "\n", inputs))
    postamble = """
    OUTFILE ${out_prefix} ${out_postfix}
    ANALYZE HETEROGENEITY
    """
    script = preamble + "\n" + processes + "\n" + postamble
    scriptFile = open("script.metal", "w")
    scriptFile.write(script)
    scriptFile.close()
    EOF
    echo "=== BEGIN metalcast.py ==="
    cat metalcast.py
    echo "=== END metalcast.py ==="
    python metalcast.py
    echo "=== BEGIN script.metal ==="
    cat script.metal
    echo "=== END script.metal ==="
    /metal script.metal
  }
  output {
    File out = glob(out_prefix + "*" + out_postfix)[0]
  }
}