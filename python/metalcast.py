scriptFile = open("script.metal", "w")
script = """
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
    """
scriptFile = open("script.metal", "w")
scriptFile.write(script)
scriptFile.close()
