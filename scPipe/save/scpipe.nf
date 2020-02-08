#!usr/bin/env nextflow

// Params
params.PlatForm = 'DropSeq'
params.seqdir = '/mnt/data/tst/Pipetool/data/Dropseq_data'
// Output file
params.outdir = './scPipe_output'
// Process Parameters
params.star_index = '/mnt/data/Static_file/dropseq_metadata/STAR_index'
params.gtf = '/mnt/data/Static_file/dropseq_metadata/mm10.gtf'
params.threads = '20'

scriptR="${workflow.projectDir}/scpipe.R"



def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run scpipe.nf --seqdir [ ] --star_index [ ]  --gtf [ ]  --PlatForm [ ] --outdir [ ] --threads [ ]

    Mandatory Arguments:

      --seqdir                      The input data directory

      --star_index                  The rapmap index directory

      --gtf                         The path of file.gtf

      --PlatForm                    Data source protocol,such as DropSeq, 10X...default is DropSeq

    Optional Arguments:
      --outdir                      The output directory where the results will be saved, default is ./scPipe_output

      --threads                     The number of threads,default is 4

    """.stripIndent()
}

params.help = null

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// Data Channel
read_pairs = Channel
	.fromPath("$params.seqdir/*_{1,2}.fastq.gz")
	.ifEmpty{error "can not find any matching: $params.seqdir"}
	.map{it->
	m=it.getName()=~/(.+)\_\d+.fastq.gz/
	tuple(m[0][1],it)
	 }
	.groupTuple(size:2, sort: true)
	.map{sample,reads->tuple(sample,reads[0],reads[1])}

process Trim_Barcode{
  cache "deep";tag "${name}.step1"
  input:
  set name, read1, read2 from read_pairs
  output:
  set name, file("trimmed.fastq") into fqTrimmed

  """
  Rscript ${scriptR} Trim_Barcode ${params.PlatForm} ${read2} ${read1} \$PWD/trimmed.fastq
  """
}

(fqTrimmed1, fqTrimmed2)=fqTrimmed.into(2)

process STAR {
  cache "deep";tag "${name}.step2"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from fqTrimmed1
  
output:
  file "${name}.Aligned.sortedByCoord.out.bam" into bam

  """
  STAR \\
  --runMode alignReads \\
  --genomeDir $params.star_index \\
  --outFileNamePrefix ${name}. \\
  --readFilesIn $x \\
  --outSAMtype BAM SortedByCoordinate \\
  --outFilterMultimapNmax 1 \\
  --runThreadN $params.threads \\
  --limitIObufferSize 1500000000 --limitOutSJcollapsed 10000000
  """
}

process Detect_Barcode{
  cache "deep";tag "${name}.step3"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'
  
  input:
  set name, x from fqTrimmed2

  output:
  set val("${name}"), file("${name}.sample_index.csv") into sampleIndex

  """
  Rscript ${scriptR} Detect_Barcode ${params.PlatForm} $x \$PWD/${name}.sample_index.csv
  """
}

process Count_Aligned_Bam{
  cache "deep";tag "${name}.step4"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from sampleIndex
  file y from bam

  output:
  set file("${name}.gene_count.csv"),file("./${name}.stat"),file("${name}.aligned.mapped.bam")

  """
  Rscript ${scriptR} Count_Aligned_Bam ${params.PlatForm} ${params.gtf} $x $y \$PWD/${name}.aligned.mapped.bam
  mv gene_count.csv ${name}.gene_count.csv
  mv ./stat ./${name}.stat
  """
}
