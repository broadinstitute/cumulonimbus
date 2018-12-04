inputs = ["yo", "blub"]
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
