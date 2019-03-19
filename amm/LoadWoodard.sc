import $ivy.`org.broadinstitute::woodard-cromiam:0.1.2`
import woodard.woosh.WoodardShell.terraProd._
import better.files.File
import ammonite.ops._
implicit def os2better(file: os.Path): File = File(file.toNIO)
