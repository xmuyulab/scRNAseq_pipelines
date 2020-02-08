#!usr/bin/env  nextflow

//Configuration
params.I = '/mnt/data/tst/Pipetool/data/10x_data'
params.O = './CellRanger_output'
params.r = '/mnt/data/Static_file/10x_cellRanger/refdata-cellranger-GRCh38-3.0.0'
params.set_cell_number = '100'
params.sample = ''

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run cellranger.nf --I [ ] --O [ ]  --r [ ]  --sample [ ] --set_cell_number [ ]

    Mandatory arguments:
      --I                      The input FASTQ data directory
      --O                      The output directory where the results will be saved, default is  ./CellRanger_output
      --r                      Path of folder containing 10x-compatible reference.
      --sample                 Prefix of the filenames of FASTQs to select.
      --set_cell_number        Expected number of recovered cells.

    """.stripIndent()
}
params.help = null
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
//params.help = null

//The main program
process cellRanger {

  cache "deep";tag "step1"
  publishDir path: "${params.O}", mode: 'copy'

  output:
  file  "*"

  """
  cellranger count \\
  --id=results \\
  --transcriptome=$params.r \\
  --fastqs=$params.I \\
  --sample=$params.sample \\
  --expect-cells=$params.set_cell_number
  """
}

// -----------------------------------------------------------------------------
workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}
