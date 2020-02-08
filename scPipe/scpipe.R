args<-commandArgs(T);
Module<-args[1];
PlatForm<-args[2];

library(stringr)
library(scPipe)

structure_DropSeq = list(bs1 = -1,bl1 = 0, bs2 = 0, bl2 = 12, us = 12, ul = 8)
structure_10X = list(bs1 = -1,bl1 = 0, bs2 = 0, bl2 = 16, us = 16, ul = 12)
read_structure=structure_DropSeq
if(PlatForm == "DropSeq"){read_structure=structure_DropSeq}
if(PlatForm == "10X"){read_structure=structure_10X}

if(Module == "Trim_Barcode"){
  fq_R1<-args[3];
  fq_R2<-args[4];
  combined_fastq<-file.path(args[5])

  sc_trim_barcode(outfq = combined_fastq, r1 = fq_R1, r2 = fq_R2, read_structure = read_structure)
}

if(Module == "Detect_Barcode"){
  combined_fastq<-args[3];
  barcode_anno<-args[4];

  CMD=paste("wc -l",combined_fastq)
  readscount=read.table(text = system(CMD,intern=T))[1]/4
  readscount=as.numeric(readscount)
  sc_detect_bc(
  infq = combined_fastq,
  outcsv = barcode_anno,
  bc_len = read_structure$bl2,
  max_reads= readscount,
  min_count = min(100,readscount/1000)
  )
}

if(Module == "Count_Aligned_Bam"){
  exon_anno<-args[3]
  barcode_anno<-args[4];
  aligned_bam<-file.path(args[5])
  mapped_bam<-file.path(args[6])

  sc_count_aligned_bam(
  inbam = aligned_bam,
  outbam = mapped_bam,
  annofn = exon_anno,
  bc_len = read_structure$bl2,
  UMI_len = read_structure$ul,
  outdir = ".",
  bc_anno = barcode_anno
   )
}
