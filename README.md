# PIPELINE
Comparison of High-Throughput Single-Cell RNA Sequencing Data Processing Pipelines

# Install Nextflow
https://www.nextflow.io/docs/latest/getstarted.html#requirements

## Requirements

1、Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 8 (or later, up to 11) to be installed.

`$ sudo apt-get install openjdk-8-jdk`

2、`$ mkdir Nextflow && cd Nextflow`

3、`$ wget -qO- https://get.nextflow.io | bash`

It will create the nextflow main executable file in the current directory.

4、move the nextflow file to a directory accessible by your $PATH variable (this is only required to avoid remembering and typing the full path to nextflow each time you need to run it).

5、Test Nextflow

`$ nextflow run hello`

Output the following to indicate that the installation was successful

```
N E X T F L O W  ~  version 19.04.1
Pulling nextflow-io/hello ...
downloaded from https://github.com/nextflow-io/hello.git
Launching `nextflow-io/hello` [nice_wozniak] - revision: a9012339ce [master]
[warm up] executor > local
executor >  local (4)
[bc/a32c37] process > sayHello [100%] 4 of 4 ✔
Ciao world!

Hola world!

Hello world!

Bonjour world!
```
# Install Anaconda
Anaconda 2019.10 for Linux Installer
Python 3.7 version

https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh

Make sure conda is installed and updated:

```
$ conda --version
$ conda update conda
```

### Tip:

If this is your first time working with Conda, you may need to edit your configuration paths to ensure Anaconda is invoked when calling `conda`

# Quickstart

## Clone Repository

`$ git clone git@github.com:xmuyulab/PIPELINE.git`

## Install Dependencies

Configuring the software environment.
Provides an execution environment for seven pieces of software, which can be installed as needed

```
$ conda env create -f PIPELINE/envs/UMI-tools.yml
$ conda env create -f PIPELINE/envs/umis.yml
$ conda env create -f PIPELINE/envs/dropEst.yml
$ conda env create -f PIPELINE/envs/zUMIs.yml
$ conda env create -f PIPELINE/envs/scPipe.yml
$ conda env create -f PIPELINE/envs/Dropseq.yml
$ conda env create -f PIPELINE/envs/cellranger.yml
```

# Drop-seq_tools

## preprocessing

### 1、Download Drop-seq_tools software

`$ wget -c https://github.com/broadinstitute/Drop-seq/releases/download/v2.3.0/Drop-seq_tools-2.3.0.zip`

`$ gunzip Drop-seq_tools-2.3.0.zip`

`$ cd Drop-seq_tools-2.3.0`
 
### 2、Modify the temporary directory of all executables, all the positions are 10 rows from the bottom. 

Attention please! 

Except file of `Drop-seq_alignment.sh`. Replace the following `/mnt/data/lmy/TMP_Dropseq` with a temporary large space directory. My working directory is `broad_tmpdir_root=/mnt/data/lmy/TMP_Dropseq`.

Please replace `/mnt/data/lmy/TMP_Dropseq`to `/path/to/tmp`，you need to create this directory

`$ mkdir -p /path/to/tmp`

### 3、Modify `Drop-seq_alignment.sh`, Add temporary directory `-Djava.io.tmpdir=/mnt/data/lmy`, please modify `/mnt/data/lmy` to `/path/to/javatmp`. And modify STAR argument. Position is located in the follow.

```
# Stage 2: alignment
$echo_prefix java -Xmx4g -Djava.io.tmpdir=/mnt/data/lmy -jar ${picard_jar} SamToFastq INPUT=${tmpdir}/unaligned_mc_tagged_polyA_filtered.bam \
  FASTQ=$tmpdir/unaligned_mc_tagged_polyA_filtered.fastq
files_to_delete="$files_to_delete $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq"

$echo_prefix $star_executable --genomeDir ${genomedir} --outFileNamePrefix ${tmpdir}/star. \
  --readFilesIn $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq --runThreadN 4 --limitIObufferSize 1500000000 --limitOutSJcollapsed 10000000
files_to_delete="$files_to_delete ${aligned_sam}"

# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
$echo_prefix java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4g -Djava.io.tmpdir=/mnt/data/lmy -jar ${picard_jar} \
  SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR=${tmpdir}
files_to_delete="$files_to_delete ${aligned_sorted_bam}"

# Stage 4: merge and tag aligned reads
$echo_prefix java -Xmx4g -Djava.io.tmpdir=/mnt/data/lmy -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${tagged_unmapped_bam} \
  ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false CLIP_ADAPTERS=false \
  TMP_DIR=${tmpdir} OUTPUT=$tmpdir/merged.bam
files_to_delete="$files_to_delete $tmpdir/merged.bam"
```
### Execute the script

4、`$ cd PIPELINE/Drop-seq_tools`

5、`--help` parameters help the user with this script.

`$ nextflow run drop_seqtools.nf --help`

6、`$ nextflow run drop_seqtools.nf --I [ ]  --gtf [ ]  --d [ ] --g [ ]  --set_cell_number [ ] --fasta [ ]  --tmpdir [ ]  --O [ ]  --javatmp [ ]`

# umis

1、`$ cd PIPELINE/umis`

2、`--help` parameters help the user with this script.

 `$ nextflow run umis.nf --help`

3、Fill in the parameters as follows

`$ nextflow run umis.nf --seqdir [ ] --rapmap_index [ ]  --set_cell_number [ ] --json [ ] --outdir [ ] --threads [ ]`

# UMI-tools

1、`$ cd PIPELINE/UMI-tools`

2、`--help` parameters help the user with this script.

`$ nextflow run umitools.nf --help`

3、`$ nextflow run umitools.nf --seqdir [ ]  --gtf [ ] --star_index [ ] --bc_pattern [ ] --set_cell_number [ ] --outdir [ ] --threads [ ]`

# dropEst

1、`$ cd PIPELINE/dropEst`

2、`--help` parameters help the user with this script.

`$ nextflow run dropEst.nf --help`

3、`$ nextflow run dropEst.nf --seqdir [ ] --star_index [ ]  --gtf [ ] --xml [ ] --outdir [ ] --threads [ ]`

# zUMIs

1、`$ cd PIPELINE/zUMIS`

2、`--help` parameters help the user with this script.

`$ nextflow run zumis.nf --help`

3、`$ nextflow run zumis.nf --shell [ ] --yaml [ ]  --outdir [ ]`

# scPipe
1、`$ cd PIPELINE/scPipe`

2、`--help` parameters help the user with this script.

`$ nextflow run scpipe.nf --help`

3、`$ nextflow run scpipe.nf --seqdir [ ] --star_index [ ]  --gtf [ ]  --PlatForm [ ] --outdir [ ] --threads [ ]`

# cellranger

1、Download cellranger

https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/

`$ wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1578322576&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU3ODMyMjU3Nn19fV19&Signature=DWyftPc5rsvMJwA3dCcHlF1ZRUgBpev3x654ICzHWke4pruaXg1D7DGcQds-Tf2-jihIWf3H3GSr9EPs-3K5dkBg9GDH61K~o5y7VVkPNEhnD~-o6NG0zNgJu-NNIh7aQE6pnBcRyMjQOVGu2~~dEt0gwdxW~GOTZe3OkJrLf8IRuwL0oEcHf9zEAxsoK8t8n9A~iS8-khzQnMtUJr0aL2Ah37mCdjRl-syfpQR~zR2InTSJkMN5TicBgBQCrfapQ6Y2q6-GAG8bHQM-aCY73d0OjdKwfKg2WFDt6hCJ5bpxCwKzwcGjaXy-pQO5YOwkGDyIX7uOv~bCcWJ7wm-ZwQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"`

2、 Decompression
`$ tar -xvf cellranger-3.0.2.tar`

3、Add environment variables

`vi ~/.bashrc`

Append 

`export PATH="/path/to/cellranger-3.0.2:$PATH"`

Save exit edit

`source ~/.bashrc`

### Execute the script

4、`$ cd PIPELINE/CellRanger`

5、`--help` parameters help the user with this script.

`$ nextflow run cellranger.nf --help`

6、`$ nextflow run cellranger.nf --I [ ] --O [ ]  --r [ ]  --set_cell_number [ ]`

## Cite


