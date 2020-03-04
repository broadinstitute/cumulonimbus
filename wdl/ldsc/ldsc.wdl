version 1.0

workflow ldsc {

}

task ldsc_task {
  command <<<
    set -e
    python3 --version
    python make_annot.py \
    		--gene-set-file GTEx_Cortex.GeneSet \
    		--gene-coord-file ENSG_coord.txt \
    		--windowsize 100000 \
    		--bimfile 1000G.EUR.QC.22.bim \
    		--annot-file GTEx_Cortex.annot.gz
    python ldsc.py\
    		--l2\
    		--bfile 1000G.EUR.QC.22\
    		--ld-wind-cm 1\
    		--annot Brain_DPC_H3K27ac.annot.gz\
    		--thin-annot
    		--out Brain_DPC_H3K27ac\
    		--print-snps hm.22.snp
    python ldsc.py
    	--h2 BMI.sumstats.gz\
    	--ref-ld-chr baseline.\
    	--w-ld-chr weights.\
    	--overlap-annot\
    	--frqfile-chr 1000G.mac5eur.\
    	--out BMI_baseline
  >>>
}