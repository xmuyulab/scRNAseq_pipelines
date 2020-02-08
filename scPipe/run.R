library(stringr)
library(scPipe)

fq_R1 = "/mnt/md0/DataProcess/SingleCellRNA/Pipetool-master/data/Dropseq_data/0.01m_2.fq.gz"
fq_R2 = "/mnt/md0/DataProcess/SingleCellRNA/Pipetool-master/data/Dropseq_data/0.01m_1.fq.gz"
out_dir = "output"

combined_fastq = file.path(out_dir, "trimmed.fastq")
aligned_bam = file.path(out_dir, "Aligned.sortedByCoord.out.bam")
mapped_bam = file.path(out_dir, "aligned.mapped.bam")


read_structure = list(
	bs1 = -1,
	bl1 = 0,
	bs2 = 0,
	bl2 = 12,
	us = 12,
	ul = 8
)


 sc_trim_barcode(
outfq = combined_fastq,
r1 = fq_R1,
r2 = fq_R2,
read_structure = read_structure
)
to_call = function(x) {
x %>% str_split(pattern = "\n") %>%
 unlist() %>% str_trim() %>% paste(collapse = " ")
}
make_call = function(x) x %>% to_call() %>% system()


 make_call(
     "STAR
         --runMode alignReads
         --runThreadN 4
         --genomeDir /mnt/md0/DataProcess/SingleCellRNA/dropseq_metadata/STAR_index
	 --readFilesIn output/trimmed.fastq
         --outSAMtype BAM SortedByCoordinate
         --outFileNamePrefix output/
         --outFilterMultimapNmax 1"
 )


barcode_anno = "output/sample_index.csv"
exon_anno="/mnt/md0/DataProcess/SingleCellRNA/dropseq_metadata/mm10.gtf"

sc_detect_bc(
	infq = combined_fastq,
	outcsv = barcode_anno,
	bc_len = read_structure$bl2,
	max_reads= 10000,
	min_count = 100
)


sc_count_aligned_bam(
inbam = aligned_bam,
outbam = mapped_bam,
annofn = exon_anno,
bc_len = read_structure$bl2,
UMI_len = read_structure$ul,
outdir = out_dir,
bc_anno = barcode_anno
 )
