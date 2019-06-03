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
if [ ! -d "genome_ind" ] 
 then  ## Genome indexes are generated if necessary 
    ./STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir /gpfs/home/juagarcia/genome_ind \
    --genomeFastaFiles /gpfs/home/juagarcia/Nitab-v4.5_genome_Scf_Edwards2017.fasta \
    --sjdbOverhang 142 \
    --sjdbGTFfile /gpfs/home/juagarcia/Galaxy64-[gffread_on_data_51__gtf].gtf \
    --genomeChrBinNbits 12 
fi   
#./STAR --runThreadN 16 --runMode alignReads --outFilterMatchNmin 16 --outSAMtype BAM Unsorted SortedByCoordinate --genomeDir /gpfs/home/juagarcia/genome_ind --outFileNamePrefix sample --readFilesIn $1  $2


## In featurecOunts, 
./featureCounts -p -C -O -t exon -g gene_id -F SAF -Q 32 -T 16 -a /gpfs/home/jumagari/Galaxy51-[Nitab-v4.5_gene_models_Scf_Edwards2017.gff].gff3  -o sample_count sampleAligned.out.sam
# | cut -f1,7- | sed 1d > $GENEMX
