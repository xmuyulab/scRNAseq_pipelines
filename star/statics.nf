#!usr/bin/env  nextflow
//results
params.outdir = "${workflow.projectDir}/genome"
params.fasta = '/mnt/data/Static_file/dropseq_metadata/mm10.fasta'
params.gtf = '/mnt/data/Static_file/dropseq_metadata/mm10.gtf'
params.tools = '/mnt/md0/Drop-seq_tools-2.3.0'

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run statics.nf --fasta [ ]  --gtf [ ]   --tools [ ] --outdir [ ]

    Mandatory arguments:
      --fasta                       /path/to/genome.fasta

      --gtf                         /path/to/genome.gtf
      --tools                       Directory containing Drop-seq executables

    Optional Arguments:
      --outdir                      The output directory where the results will be saved. Default: current directory ./genome

    """.stripIndent()
}
params.help = null

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

process CreateSeqDict {

  cache "deep";tag "step1"
  publishDir path: "${params.outdir}", mode: 'copy'

  output:
  file  "genome.dict" into refFlat

  """
  gatk CreateSequenceDictionary -R $params.fasta -O genome.dict
  """
}

(RGtf, RefFlat, CIFiles1)=refFlat.into(3)

process ConvertToRefFlat {

  cache "deep";tag "step2"
  publishDir path: "${params.outdir}", mode: 'copy'

  input:
  file x from RefFlat

  output:
  file  "genome.refFlat" 

  """
  $params.tools/ConvertToRefFlat ANNOTATIONS_FILE=$params.gtf SEQUENCE_DICTIONARY=$x OUTPUT=genome.refFlat
  """
}

process ReduceGtf {

  cache "deep";tag "step3"
  publishDir path: "${params.outdir}", mode: 'copy'

  input:
  file x from RGtf

  output:
  file  "genome.reduced.gtf" into CIFiles2

  """
  $params.tools/ReduceGtf SEQUENCE_DICTIONARY=$x GTF=$params.gtf OUTPUT=genome.reduced.gtf
  """
}

process CInFiles {

  cache "deep";tag "step4"
  publishDir path: "${params.outdir}", mode: 'copy'

  input:
  file x from CIFiles1
  file y from CIFiles2

  output:
  file  "*" 

  """
  $params.tools/CreateIntervalsFiles SEQUENCE_DICTIONARY=$x REDUCED_GTF=$y PREFIX=genome OUTPUT=./ MT_SEQUENCE=MT
  """
}


// -----------------------------------------------------------------------------
workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}


