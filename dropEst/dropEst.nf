#!usr/bin/env  nextflow
// nextflow run dropEst.nf -with-report dropEst.html -with-dag dropEst.png
/* Process Parameters */

// Input
params.seqdir = '/mnt/data/tst/Pipetool/data/Dropseq_data'
// Output
params.outdir = './dropEst_output'

params.star_index = '/mnt/data/Static_file/dropseq_metadata/STAR_index'
params.gtf = '/mnt/data/Static_file/dropseq_metadata/mm10.gtf'
params.threads = 10
params.xml = "/mnt/data/tst/1225/pipelines/dropEst/configs/drop_seq.xml"
//params.xml = null


def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run dropEst.nf --seqdir [ ] --star_index [ ]  --gtf [ ] --xml [ ] --outdir [ ] --threads [ ]

    Mandatory Arguments:

      --seqdir                      The input data directory

      --star_index                  The rapmap index directory

      --gtf                         The path of file.gtf

      --xml                         config filename, xml file with droptag parameters.

    Optional Arguments:
      --outdir                      The output directory where the results will be saved, default is ./dropEst_output

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

// Config Channel
//xml = Channel.fromPath("${workflow.projectDir}/configs/drop_seq.xml")
xml = Channel.fromPath("${params.xml}")

process Droptag {

  cache "deep";tag "${name}.step1"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:// nextflow里，若直接接受多个channel的输入，则pipeline会将这些turple中的element逐个排列，进入Process
	// 因此，若实现多个Channel的所有组合作为输入产生下一步结果，就必须要用combine的option进行了。
	// 另外，若多个组合，需要进行诸如配对的组合时，务必采用combine这类的option，防止出现意外
  set name, read1, read2, x from read_pairs.combine(xml)

  output:
  set val("${name}"), file("${name}_2.fastq.gz.tagged.fastq.gz") into to_star
  set file("${name}_2.fastq.gz.tagged.params.gz"), val("${x}") into to_next
  file "*"

  """
  droptag -c $x -r 0 -p $params.threads -S -s $read1 $read2
  """
}

process STAR {

  cache "deep";tag "${name}.step2"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from to_star

  output:
  set val("${name}"), file ("${name}.Aligned.out.bam") into to_dropest
  file ("${name}.*")

  """
  STAR \\
  --genomeDir $params.star_index \\
  --readFilesIn $x \\
  --outFileNamePrefix ${name}. \\
  --outSAMtype BAM Unsorted \\
  --readFilesCommand zcat \\
  --runThreadN $params.threads \\
  --limitIObufferSize 1500000000 --limitOutSJcollapsed 10000000
  """
}


process Dropest {

  cache "deep";tag "${name}.step3"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set x, z from to_next
  set name, y from to_dropest
  output:
  // set val("${name}"), file ("${name}.cell.counts.rds") into follow
  file ("${name}.*")

  """
  dropest -w -m -r $x -g $params.gtf -c $z $y
  ls *|xargs -i mv {} ${name}.{}
  """
}
// -----------------------------------------------------------------------------

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}

