#!usr/bin/env nextflow
// nextflow run umis.nf -with-report umis.html -with-dag umis.png
// Input
params.seqdir = '/mnt/data/tst/Pipetool/data/Dropseq_data'
// Output
params.outdir = './umis_output'
// Process Parameters
params.rapmap_index = '/mnt/data/txw/umis/mm10_data_runtime/assistant_files/ref_index'
params.threads = '30'
params.set_cell_number = '100'

params.json = "/mnt/data/tst/1225/pipelines/umis/configs/Dropseq_transform.json"
//params.json = null

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run umis.nf --seqdir [ ] --rapmap_index [ ]  --set_cell_number [ ] --json [ ] --outdir [ ] --threads [ ]

    Mandatory Arguments:

      --seqdir                      The input data directory

      --rapmap_index                The rapmapindex directory

      --set_cell_number             The number of cells
      
      --json                        The command umis fastqtransform is for transforming a (pair of) read(s) to this format based on a transform file. 
                                    The transform file is a json file which has a Python flavored regular expression for each read, 
                                    made to extract the necessary components of the reads. 

    Optional Arguments:
      --outdir                      The output directory where the results will be saved, default is ./umis_output

      --threads                     The number of threads,default is 4

    """.stripIndent()
}

params.help = null

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


read_pairs = Channel
        .fromPath("$params.seqdir/*_{1,2}.fastq.gz")
        .ifEmpty{error "can not find any matching: $params.seqdir"}
        .map{it->
        m=it.getName()=~/(.+)\_\d+.fastq.gz/
        tuple(m[0][1],it)
         }
        .groupTuple(size:2, sort: true)
        .map{sample,reads->tuple(sample,reads[0],reads[1])}

//json = Channel.fromPath("${workflow.projectDir}/configs/Dropseq_transform.json")
json = Channel.fromPath("${params.json}")

process Fastqtransform {

  cache "deep";tag "${name}.step1"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:// nextflow里，若直接接受多个channel的输入，则pipeline会将这些turple中的element逐个排列，进入Process
	// 因此，若实现多个Channel的所有组合作为输入产生下一步结果，就必须要用combine的option进行了。
	// 另外，若多个组合，需要进行诸如配对的组合时，务必采用combine这类的option，防止出现意外
  set name, read1, read2 , x from read_pairs.combine(json)

  output:
  set val("${name}"),file ("${name}.formatted.fq") into to_histogram
  set val("${name}"),file ("${name}.formatted.fq") into to_quasimap

  """
  umis fastqtransform $x $read1 $read2 > ${name}.formatted.fq
  """
}

process cb_histogram {

  cache "deep";tag "${name}.step2"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from to_histogram

  output:
  set val("${name}"), file ("${name}.cb-histogram.txt") into to_select

  """
  umis cb_histogram $x > ${name}.cb-histogram.txt
  """
}

process rapmap {

  cache "deep";tag "${name}.step3"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from to_quasimap

  output:
  set val("${name}"), file ("${name}.final.sam") into to_bam

  """
  rapmap quasimap -i $params.rapmap_index -r $x -t $params.threads -o ${name}.final.sam
  """
}

process Samtools {

  cache "deep";tag "${name}.step4"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from to_bam

  output:
  set val("${name}"), file ("${name}.final.bam") into samtoolsort

  """
  samtools view -bS $x > ${name}.final.bam
  """
}

process SamtoolSort {

  cache "deep";tag "${name}.step5"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from samtoolsort

  output:
  set val("${name}"), file ("${name}.sorted.bam") into countone

  """
  samtools sort $x -o ${name}.sorted.bam
  """
}

process SelectHistogram {

  cache "deep";tag "${name}.step6"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from to_select

  output:
  set val("${name}"), file ("${name}.selected-cb-histogram.txt") into tagcount

  """
  head -n $params.set_cell_number $x > ${name}.selected-cb-histogram.txt
  """
}

process tagcount {

  cache "deep";tag "${name}.step7"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from countone
  set name, y from tagcount

  output:
  set val("${name}"), file ("${name}.final_result.txt") into follow

  """
  samtools index $x
  umis fasttagcount --cb_histogram $y $x ${name}.final_result.txt
  """
}

// -----------------------------------------------------------------------------

//workflow.onComplete {
//  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
//}

if(workflow.success) {

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
