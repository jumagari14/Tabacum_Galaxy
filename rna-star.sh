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

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib

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
list=$(ls subdata/*.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

mkdir -p -m 755 STAR_Align
mkdir -p -m 755 counts
for I in $list
do
    ./STAR \
    --runThreadN 16 \
    --runMode alignReads \
    --outFilterMatchNmin 16 \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --genomeDir /gpfs/home/juagarcia/genome_ind \
    --outFileNamePrefix STAR_Align/"$I" \
    --readFilesIn subdata/"$I"_1.fq.gz  subdata/"$I"_2.fq.gz

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

cd counts ; touch Matrix_data 

count=0

for J in $list 
do 
    count=$((count+1))
    if [$count -eq 1 ]
    then 
        sed '1d' $J > tmp ; mv tmp $J 
        cat $J | cut -d $'\t' -f 1,7  > count_ind
        paste Matrix_data count_ind > temp2 && mv temp2 Matrix_data.tabular 
        rm -f temp2
    else
        sed '1d' $J > tmp ; mv tmp $J 
        cat $J | rev | cut -d $'\t' -f 1 | rev > count_ind
        paste Matrix_data count_ind > temp2 && mv temp2 Matrix_data.tabular 
        rm -f temp2 
    fi 
done 

rm -f count_ind 
