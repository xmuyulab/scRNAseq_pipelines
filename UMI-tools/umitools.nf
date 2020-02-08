#!usr/bin/env  nextflow
// nextflow run umitools.nf -with-report umitools.html -with-dag umitools.png 
// Input file
params.seqdir = '/mnt/data/tst/Pipetool/data/Dropseq_data'
// Output file
params.outdir = './UMItools_output'
// Process Parameters
params.star_index = '/mnt/data/Static_file/dropseq_metadata/STAR_index'
params.gtf = '/mnt/data/Static_file/dropseq_metadata/mm10.gtf'
params.bc_pattern = 'CCCCCCCCCCCCNNNNNNNN'
params.set_cell_number = '100'
params.threads = '20'


def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    \$ nextflow run umitools.nf --seqdir [ ]  --gtf [ ] --star_index [ ] --bc_pattern [ ] --set_cell_number [ ] --outdir [ ] --threads [ ]

    Mandatory arguments:
      --seqdir                      The input data directory

      --gtf                         The path of file.gtf

      --star_index                  The STAR index directory

      --bc_pattern                  CCCCCCCCCCCCNNNNNNNN
                                    Use C characters to show where CellBarcode bases are and N characters to show were UMI bases are.
                                    Thus, in the above we have 12 Cs followed by 8 Ns to denote 
                                    that the first 12 bases of the read are CellBarcode bases and the second 8 are UMI bases.

      --set_cell_number             The number of cells

    Optional Arguments:
      --outdir                      The output directory where the results will be saved, default is ./UMItools_output

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

(to_first, to_two)=read_pairs.into(2)

process Whitelist {

  cache "deep";tag "${name}.step1"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, read1, read2 from to_first

  output:
  file "${name}.whitelist.txt" into to_extract

  """
  umi_tools whitelist \\
  --stdin $read1 \\
  --bc-pattern=$params.bc_pattern \\
  --set-cell-number=$params.set_cell_number \\
  --log2stderr > ${name}.whitelist.txt
  """
}

process Extract {

  cache "deep";tag "${name}.step2"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, read1, read2 from to_two
  file x from to_extract

  output:
  set val("${name}") ,file("${name}.2_extracted.fastq.gz") into to_star
  file "${name}.1_extracted.fastq.gz"

  """
  umi_tools extract \\
  --bc-pattern=$params.bc_pattern \\
  --stdin $read1 \\
  --stdout ${name}.1_extracted.fastq.gz \\
  --read2-in $read2 \\
  --read2-out=${name}.2_extracted.fastq.gz \\
  --filter-cell-barcode \\
  --whitelist=$x
  """
}

process STAR {

  cache "deep";tag "${name}.step3"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from to_star

  output:
  set val("${name}") ,file("${name}.Aligned.sortedByCoord.out.bam") into featureCount
  file ("${name}.*")

  """
  STAR \\
  --genomeDir $params.star_index \\
  --readFilesIn $x \\
  --outFileNamePrefix ${name}. \\
  --outSAMtype BAM SortedByCoordinate \\
  --readFilesCommand zcat \\
  --outFilterMultimapNmax 1 \\
  --runThreadN $params.threads \\
  --limitIObufferSize 1500000000 --limitOutSJcollapsed 10000000
  """
}

process FeatureCount {

  cache "deep";tag "${name}.step4"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from featureCount

  output:
  file "${name}.gene_assigned*"
  set val("${name}"), file("${name}.Aligned.sortedByCoord.out.bam.featureCounts.bam") into samtoolsort

  """
  featureCounts -a $params.gtf -o ${name}.gene_assigned -R BAM $x -T $params.threads
  """
}

process Samtools {

  cache "deep";tag "${name}.step5"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x  from samtoolsort

  output:
  set val("${name}"),file("${name}.assigned_sorted.bam") into count

  """
  samtools sort $x -o ${name}.assigned_sorted.bam
  """
}


process Count {

  cache "deep";tag "${name}.step6"
  publishDir path: "${params.outdir}/${name}", mode: 'copy'

  input:
  set name, x from count

  output:
  file "${name}.counts.tsv.gz"

  """
  samtools index $x
  umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I $x -S ${name}.counts.tsv.gz
  """
}
// -----------------------------------------------------------------------------

workflow.onComplete {
  println (workflow.success ? "Success: Pipeline Completed!" : "Error: Something went wrong.")
}
