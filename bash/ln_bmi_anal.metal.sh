#!/usr/bin/env bash

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
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_sa/t2d_sa.sl_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_sa/t2d_sa.sp_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_sa/t2d_sa.ss_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_sa/t2d_sa.ss_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_ea/t2d_ea.eh_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_ea/t2d_ea.ek_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_ea/t2d_ea.ek_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_ea/t2d_ea.es_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_ea/t2d_ea.es_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_hs/t2d_hs.ha_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_hs/t2d_hs.hs_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_hs/t2d_hs.hs_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_hs/t2d_hs.sigma_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_hs/t2d_hs.sigma_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_aa/t2d_aa.aj_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_aa/t2d_aa.aw_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_aa/t2d_aa.aa_esp.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_aa/t2d_aa.ab_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.eu_esp.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.got2d.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.lucamp.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.ua_13k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.uf_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.ug_29k.ln_bmi_anal.vassoc
PROCESS /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/t2d_eu/t2d_eu.um_13k.ln_bmi_anal.vassoc
OUTFILE /home/unix/flannick/links/t2d_exomes/amp/55k/TM6SF2/vassoc/ln_bmi_anal.metal.vassoc .tbl
ANALYZE HETEROGENEITY
EOF
/metal script.metal

