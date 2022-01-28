#!/bin/bash
#SBATCH -n 16
#SBATCH --mem 100G
#SBATCH --job-name=HiCAlign 
#SBATCH --output=align_bwa_%j.txt
#SBATCH  -p hpc5,Lewis
#SBATCH -t 2-00:00:00


export MYLOCAL=/storage/hpc/data/fbqc9/.conda/envs
module load miniconda3
source activate $MYLOCAL/rna_seq

# fastq1="/storage/htc/bdm/Frimpong/hic_ana/raw_data/DS_Pool_$1_R1.fastq"
# fastq2="/storage/htc/bdm/Frimpong/hic_ana/raw_data/DS_Pool_$1_R2.fastq" 
# ref="/storage/htc/bdm/Frimpong/ref_geneomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"
# chrom_sizes="/storage/htc/bdm/Frimpong/ref_geneomes/Mus_musculus/UCSC/mm10/mm10.chrom.sizes"

# echo $fastq1
# echo $fastq2
# echo $ref

# mkdir $1

echo "Aligning $1"
# bwa mem -SP5M -t16 $ref $fastq1 $fastq2 | samtools view -bhS - > $1/$1.bam

echo "Parsing to Pairtools"
# samtools view -h $1/$1.bam | pairtools parse -c mm10.chrom.sizes -o $1/$1_parsed.pairsam.gz

echo "Sorting Pairtools"
# pairtools sort --tmpdir=./ --nproc 16 -o $1/$1_sorted.pairsam.gz $1/$1_parsed.pairsam.gz

echo "Removing Dedups"
# pairtools dedup --output-stats $1/$1_stats_dedup.pairsam.gz --mark-dups -o $1/$1_deduped.pairsam.gz $1/$1_sorted.pairsam.gz

echo "Selecting"
# pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' -o $1/$1_filtered.pairsam.gz $1/$1_deduped.pairsam.gz

echo "Split"
# pairtools split --output-pairs $1/$1_output.pairs.gz $1/$1_filtered.pairsam.gz

echo "Finished aligning $1 at  $(date) "

echo "generating pair indicies"
# pairix $1/$1_output.pairs.gz

echo "generating cooler files"

echo "40kb 50kb 80kb 1mb"
# mkdir $1/40kb
# mkdir $1/50kb
# mkdir $1/1mb
mkdir $1/80kb
mkdir $1/100kb

# cooler cload pairix mm10.chrom.sizes:40000 $1/$1_output.pairs.gz $1/40kb/$1_40kb.cool
cooler cload pairix mm10.chrom.sizes:100000 $1/$1_output.pairs.gz $1/100kb/$1_100kb.cool
# cooler cload pairix mm10.chrom.sizes:500000 $1/$1_output.pairs.gz $1/500kb/$1_50kb.cool
# cooler cload pairix mm10.chrom.sizes:800000 $1/$1_output.pairs.gz $1/800kb/$1_80kb.cool
# cooler cload pairix mm10.chrom.sizes:1000000 $1/$1_output.pairs.gz $1/1mb/$1_1mb.cool

echo "balancing cooler files"
# cooler balance $1/40kb/$1_40kb.cool
cooler balance $1/100kb/$1_100kb.cool
# cooler balance $1/50kb/$1_50kb.cool
# cooler balance $1/80kb/$1_80kb.cool
# cooler balance $1/1mb/$1_1mb.cool

echo "balancing cooler files"
# cooler zoomify $1/40kb/$1_40kb.cool
# cooler zoomify $1/50kb/$1_50kb.cool
cooler zoomify $1/100kb/$1_100kb.cool
# cooler zoomify $1/1mb/$1_1mb.cool

echo "Finished at  $(date)"





