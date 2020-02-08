#!usr/bin/env  nextflow
// nextflow run drop_seqtools.nf -with-report report.html  -with-dag Dropseq.png

// Process Parameter
params.I = '/mnt/data/tst/Pipetool/data/Dropseq_data'
params.gtf = "/mnt/data/Static_file/dropseq_metadata/mm10.gtf"
params.d ='/mnt/md0/Drop-seq_tools-2.3.0'
params.g = "/mnt/data/Static_file/dropseq_metadata/STAR_index"
params.set_cell_number = '100'
params.fasta = "/mnt/data/Static_file/dropseq_metadata/mm10.fasta"

params.tmpdir = '\$PWD/TMP'
params.javatmp = '\$PWD/Tmp'

// Output
params.O = './Dropseq_output'

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run drop_seqtools.nf --I [ ]  --gtf [ ]  --d [ ] --g [ ]  --set_cell_number [ ] --fasta [ ]  --tmpdir [ ]  --O [ ]  --javatmp [ ]

    Mandatory arguments:
      --I                           The input data directory

      --gtf                         The path of file.gtf

      --d                           Directory containing Drop-seq executables.

      --g                           Directory of STAR genome directory

      --set_cell_number             The number of cells

      --fasta                       The path of file.fasta, ask for the .fasta suffix, can not .fa 

      --tmpdir                      Where to write temporary files. default is ./TMP
    Optional Arguments:
      --O                           The output directory where the results will be saved, default is ./Dropseq_output

      --javatmp                     JAVA temporary storage directory, default is ./Tmp

    """.stripIndent()
}
params.help = null

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//params.statics ='/mnt/md0/DataProcess/SingleCellRNA/dropseq_metadata'

read_pairs = Channel
	.fromPath("$params.I/*_{1,2}.fastq.gz")
	.ifEmpty{error "can not find any matching: $params.I"}
	.map{it->
	m=it.getName()=~/(.+)\_\d+.fastq.gz/
	tuple(m[0][1],it)
	 }
	.groupTuple(size:2, sort: true)
	.map{sample,reads->tuple(sample,reads[0],reads[1])}

process FastqToSam {

  cache "deep";tag "${name}.step1"
  publishDir path: "${params.O}/${name}", mode: 'copy'

  input:
  set name, read1, read2 from read_pairs

  output:
  set val("${name}"),file("${name}.unmapped-queryname-sorted.bam") into to_first

  """
  gatk FastqToSam \\
  --java-options "-Djava.io.tmpdir=$params.javatmp" \\
  -F1 $read1 \\
  -F2 $read2 \\
  -O ${name}.unmapped-queryname-sorted.bam \\
  -SM 0_unmapped-queryname-sorted
  """
}

process Dropalignment {

  cache "deep";tag "${name}.step2"
  publishDir path: "${params.O}/${name}", mode: 'copy'

  input:
  set name, x from to_first

  output:
  set val("${name}"),file("${name}.final.bam"),file("${name}.final.bai") into bam
  file "*"

  """
  mkdir $params.tmpdir
  ${params.d}/Drop-seq_alignment.sh \\
  -g $params.g \\
  -r $params.fasta  \\
  -t $params.tmpdir \\
  $x
  rm -rf $params.tmpdir
  ls *|xargs -i mv {} ${name}.{}
  """
}

(bam1,bam2)=bam.into(2)

process BAMTagHistogram {

  cache "deep";tag "${name}.step3"
  publishDir path: "${params.O}/${name}", mode: 'copy'

  input:
  set name,x,y from bam1

  output:
  file "${name}.cell_readcounts.txt.gz" into next

  """
  ${params.d}/BamTagHistogram \\
  I=$x \\
  O=${name}.cell_readcounts.txt.gz \\
  TAG=XC
  """
}

process DigitalExpression {

  cache "deep";tag "${name}.step4"
  publishDir path: "${params.O}/${name}", mode: 'copy'

  input:
  set name,x,y from bam2

  output:
  //file "2_gene_exon_tagged.dge.txt.gz" into follow
  file "*"

  """
  ${params.d}/DigitalExpression \\
  I=$x \\
  O=${name}.gene_exon_tagged.dge.txt.gz \\
  SUMMARY=${name}.out_gene_exon_tagged.dge.summary.txt \\
  NUM_CORE_BARCODES=$params.set_cell_number
  """
}

// -----------------------------------------------------------------------------
workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}
