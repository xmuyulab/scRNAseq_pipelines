#!usr/bin/env  nextflow

//Configuration
params.outdir = './zUMIs_output'
params.shell = '/mnt/data/lmy/zUMIs/zUMIs-master.sh'
//Shell = Channel.fromPath( '/mnt/data/lmy/zUMIs/zUMIs-master.sh' )
params.yaml = 'configs/runExample.yaml'
//Yaml = Channel.fromPath( 'runExample.yaml' )

def helpMessage() {
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    \$ nextflow run zumis.nf --shell [ ] --yaml [ ]  --outdir [ ] 

    Mandatory arguments:

      --shell                   Where is the ZUMIs-master.sh

      --yaml                    Path to the YAML config file

    Optional arguments

      --outdir                  The output directory where the results will be saved, default is  ./zUMIs_output

    """.stripIndent()
}

params.help = null

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

Shell = Channel.fromPath("${params.shell}")
Yaml = Channel.fromPath("${params.yaml}")

process zUMIs {

  cache "deep";tag "step1"
  publishDir path: "${params.outdir}", mode: 'copy'

  input:
  file x from Shell
  file y from Yaml

  output:
  file  "*"

  """
  bash $x -y $y
  """
}

// -----------------------------------------------------------------------------
workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}

