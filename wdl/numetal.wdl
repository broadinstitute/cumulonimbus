workflow numetal {
  File in00
  File in01
  File in02
  File in03
  File in04
  File in05
  File in06
  File in07
  File in08
  File in09
  File in10
  File in11
  File in12
  File in13
  File in14
  File in15
  File in16
  File in17
  File in18
  File in19
  File in20
  File in21
  File in22
  File in23
  File in24
  String out_prefix
  String out_postfix

  call metal {
    input:
      in00 = in00,
      in01 = in01,
      in02 = in02,
      in03 = in03,
      in04 = in04,
      in05 = in05,
      in06 = in06,
      in07 = in07,
      in08 = in08,
      in09 = in09,
      in10 = in10,
      in11 = in11,
      in12 = in12,
      in13 = in13,
      in14 = in14,
      in15 = in15,
      in16 = in16,
      in17 = in17,
      in18 = in18,
      in19 = in19,
      in20 = in20,
      in21 = in21,
      in22 = in22,
      in23 = in23,
      in24 = in24,
      out_prefix = out_prefix,
      out_postfix = out_postfix
  }

  output {
    File out = metal.out
  }
}

task metal {
  File in00
  File in01
  File in02
  File in03
  File in04
  File in05
  File in06
  File in07
  File in08
  File in09
  File in10
  File in11
  File in12
  File in13
  File in14
  File in15
  File in16
  File in17
  File in18
  File in19
  File in20
  File in21
  File in22
  File in23
  File in24
  String out_prefix
  String out_postfix
  runtime {
    docker: "us.gcr.io/broad-gdr-dig-storage/metal:2018-08-28"
    cpu: 1
    memory: "3 GB"
    disks: "local-disk 20 HDD"
  }
  command {
    cat << EOF > script.metal
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
    PROCESS ${in00}
    PROCESS ${in01}
    PROCESS ${in02}
    PROCESS ${in03}
    PROCESS ${in04}
    PROCESS ${in05}
    PROCESS ${in06}
    PROCESS ${in07}
    PROCESS ${in08}
    PROCESS ${in09}
    PROCESS ${in10}
    PROCESS ${in11}
    PROCESS ${in12}
    PROCESS ${in13}
    PROCESS ${in14}
    PROCESS ${in15}
    PROCESS ${in16}
    PROCESS ${in17}
    PROCESS ${in18}
    PROCESS ${in19}
    PROCESS ${in20}
    PROCESS ${in21}
    PROCESS ${in22}
    PROCESS ${in23}
    PROCESS ${in24}
    OUTFILE ${out_prefix} ${out_postfix}
    ANALYZE HETEROGENEITY
    EOF
    /metal script.metal
  }
  output {
    File out = glob(out_prefix + "*" + out_postfix)[0]
  }
}