#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J HelloWorld

# walltime (hh:mm::ss)
#SBATCH -t 00:10:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 1
#SBATCH --tasks-per-node=1

# change working directory
# SBATCH --chdir=.

# fin des directives PBS
#############################

# On terminal, paste sbatch -J STAR-Align -t 03:00:00 -N 16 --ntasks=16 --tasks-per-node=30 --chdir=. ./rna-star.sh

# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################

## Load necessary modules
module load R/3.5.0 

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib:/usr/bin/java:/usr/lib/java:/etc/java:/usr/share/java:/usr/share/man/man1/java.1.gz


module load Fiji/java-8
## java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /media/jumagari/JUANMA/Stage/Galaxy_An/subdata/TAB0.3_1.fastqsanger /media/jumagari/JUANMA/Stage/Galaxy_An/subdata/TAB0.3_2.fastqsanger forw_par.gq.gz forw_unp.fq.gz rev_pair.fq.gz rev_unp.fq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50

list=$(ls subdata/*.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

mkdir -p -m 755 trimm_data

gtf=$(find ./ -name "*.gtf")
gff=$(find ./ -name "*gff[3]*")
fasta=$(find ./ -name "*.fasta")

if [ ! -d "genome_ind" ] 
 then  ## Genome indexes are generated if necessary 
    mkdir -m 755 -p genome_ind 
    ./STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir /gpfs/home/juagarcia/genome_ind \
    --genomeFastaFiles "$fasta" \
    --sjdbOverhang 142 \
    --sjdbGTFfile "$gtf" \
    --genomeChrBinNbits 12 
fi  

mkdir -p -m 755 STAR_Align
mkdir -p -m 755 counts
for I in $list
do

    cp subdata/"$I"_1.fq.gz subdata/"$I"_1.fastqsanger 
    cp subdata/"$I"_2.fq.gz subdata/"$I"_2.fastqsanger

    java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 ./subdata/"$I"_1.fastqsanger ./subdata/"$I"_2.fastqsanger trimm_data/"$I"_1.par.fastqsanger "$I"_1.unp.fastqsanger trimm_data/"$I"_2.par.fastqsanger "$I"_2.unp.fastqsanger ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50

    rm -rf "$I"_1.unp.fastqsanger "$I"_2.unp.fastqsanger

    ./STAR \
    --runThreadN 16 \
    --runMode alignReads \
    --outFilterMatchNmin 16 \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --genomeDir /gpfs/home/juagarcia/genome_ind \
    --outFileNamePrefix STAR_Align/"$I" \
    --readFilesIn trimm_data/"$I"_1.par.fastqsanger  trimm_data/"$I"_2.par.fastqsanger

    ## In featureCounts, paired-end reads must have -p option!!! 
    ./featureCounts \
    -p \
    -t exon \
    -g gene_id \
    -F GTF \
    -Q 32 \
    -T 16 \
    -a "$gtf" \
    -o /gpfs/home/juagarcia/counts/"$I" \
    STAR_Align/"$I"Aligned.out.bam 
    # | cut -f1,7- | sed 1d > $GENEMX
done

cd counts ; touch Matrix_data.tabular 

count=0

for J in $list 
do 
    count=$((count+1))
    if (($count == 1))
    then 
        cut -d $'\t' -f 1,7 $J  > count_ind
        paste Matrix_data.tabular count_ind > temp2 && mv temp2 Matrix_data.tabular 
        rm -f temp2
    else
        cut -d $'\t' -f 7 $J  > count_ind
        paste Matrix_data.tabular count_ind > temp2 && mv temp2 Matrix_data.tabular 
        rm -f temp2 
    fi 
done 


rm -f count_ind 

## Rscript edgeR.r
## Rscript deseq2.r 