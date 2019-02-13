in_file = argv[1]
out_file = argv[2]
dat = read.table(in_file, header=T)
jpeg(out_file, width=1200, height=1000)
plot(-log(dat$P1), -log(dat$P2))
dev.off()
