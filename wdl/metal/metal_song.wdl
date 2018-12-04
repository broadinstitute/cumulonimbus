import "metal.wdl"

workflow metal_song {
  Array[File] inputs
  String out_prefix
  String out_postfix

  call metal.metal {
    input:
      inputs = inputs,
      out_prefix = out_prefix,
      out_postfix = out_postfix
  }

  output {
    File out = metal.out
  }
}