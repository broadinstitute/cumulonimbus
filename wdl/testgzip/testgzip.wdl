workflow testgzip {
 File in

 call gzip {
   input:
     in = in,
     out_name = basename(in) + ".gz"
 }

 output {
   File gzipped = gzip.out
 }
}

task gzip {
 File in
 String out_name
 command {
   gzip ${in} -c >${out_name}
 }
 runtime {
   docker: "ubuntu:16.04"
   cpu: 1
   memory: "3 GB"
   disks: "local-disk 20 HDD"
 }
 output {
   File out = out_name
 }
}

