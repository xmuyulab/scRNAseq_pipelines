#!usr/bin/env  nextflow
//results
params.outdir = "${workflow.projectDir}"
params.fasta = '/mnt/data/Static_file/ERCC/genome.fasta'
params.gtf = '/mnt/data/Static_file/ERCC/genome.gtf'
params.threads = '30'

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run star_index.nf --fasta [ ]  --gtf [ ]  --outdir [ ] --threads [ ]

    Mandatory arguments:
      --fasta                       /path/to/genome.fasta

      --gtf                         /path/to/genome.gtf

    Optional Arguments:
      --outdir                      The output directory where the results will be saved. Default: current directory ./STAR_index

      --threads                     The number of threads, default is 4

    """.stripIndent()
}
params.help = null

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

process StarIndex {

  cache "deep";tag "step1"
  publishDir path: "${params.outdir}", mode: 'copy'

  output:
  file  "*"

  """
  mkdir STAR_index
  STAR \\
  --runThreadN ${params.threads} \\
  --runMode genomeGenerate \\
  --genomeDir STAR_index \\
  --genomeFastaFiles ${params.fasta} \\
  --sjdbGTFfile ${params.gtf}
  """
}

// -----------------------------------------------------------------------------
workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}


