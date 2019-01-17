dat <- read.table("~/bottomLine/results/joined/pValuesBy1024", header=T)
plot(-log(dat$P1), -log(dat$P2))
