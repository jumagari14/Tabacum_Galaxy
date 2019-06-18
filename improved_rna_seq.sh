#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J RNA-Seq analysis

# walltime (hh:mm::ss)
#SBATCH -t 00:30:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 2
#SBATCH --tasks-per-node=16

# change working directory
# SBATCH --chdir=.

# fin des directives PBS
#############################

# On terminal, paste sbatch -J STAR-Align -t 03:00:00 --ntasks=32 --chdir=. ./improved_rna_seq.bash 

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
module load jdk1.8/8u22

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib

## java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /media/jumagari/JUANMA/Stage/Galaxy_An/subdata/TAB0.3_1.fq.gz /media/jumagari/JUANMA/Stage/Galaxy_An/subdata/TAB0.3_2.fq.gz forw_par.gq.gz forw_unp.fq.gz rev_pair.fq.gz rev_unp.fq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50

list=$(ls subdata/*.fq | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

mkdir -p -m 755 trimm_data 
mkdir -p -m 755 trimm_data/quality

gtf=$(find ./ -not -path '*/\R*' -name "*.gtf")
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

## Only to be run one time 
# if [$1=="convert"] 
# then 
#     parallel -j 3 "mv subdata/{}_1.fq.gz subdata/{}_1.fastqsanger" ::: $list 
#     parallel -j 3 "mv subdata/{}_2.fq.gz subdata/{}_2.fastqsanger" ::: $list
# fi

# Uncompress gz files 
# parallel -j 3 "gunzip {}_1.fq.gz" ::: $list

mkdir -p -m 755 STAR_Align  
mkdir -p -m 755 counts

#parallel -j 6 "FastQC/fastqc subdata/{}_1.fq subdata/{}_2.fq --outdir=trimm_data/quality" ::: $list

parallel -j 6 "java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 ./subdata/{}_1.fq ./subdata/{}_2.fq trimm_data/{}_1.par.fq {}_1.unp.fq trimm_data/{}_2.par.fq {}_2.unp.fq ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50" ::: $list

parallel -j 6 "rm -rf {}_1.unp.fq {}_2.unp.fq" ::: $list
# for I in $list
# do
# ./STAR \
#     --runThreadN 16 \
#     --runMode alignReads \
#     --outFilterMatchNmin 16 \
#     --outSAMtype BAM Unsorted SortedByCoordinate \
#     --genomeDir /gpfs/home/juagarcia/genome_ind \
#     --outFileNamePrefix STAR_Align/"$I" \
#     --readFilesIn trimm_data/"$I"_1.par.fq  trimm_data/"$I"_2.par.fq
# done
# ## In featureCounts, paired-end reads must have -p option!!! 
# cd STAR_Align
# parallel -j 3 "../featureCounts -p -t exon -g gene_id -F GTF -Q 32 -T 16 -a ../"$gtf" -o /gpfs/home/juagarcia/counts/{.}.txt {}" ::: *Aligned.out.bam 
# # | cut -f1,7- | sed 1d > $GENEMX

# cd ~/counts ; touch Matrix_data.tabular 

# count=0

# for J in $list 
# do 
#     count=$((count+1))
#     if (($count == 1))
#     then 
#         sed '1d' "$J"Aligned.out.txt > tmp ; mv tmp "$J"Aligned.out.txt 
#         cut -d $'\t' -f 1,7 "$J"Aligned.out.txt  > count_ind
#         paste Matrix_data.tabular count_ind > temp2 && mv temp2 Matrix_data.tabular 
#         rm -f temp2
#     else
#         sed '1d' "$J"Aligned.out.txt > tmp ; mv tmp "$J"Aligned.out.txt 
#         cut -d $'\t' -f 7 "$J"Aligned.out.txt  > count_ind
#         paste Matrix_data.tabular count_ind > temp2 && mv temp2 Matrix_data.tabular 
#         rm -f temp2 
#     fi 
# done 

# cd ~/; rm -f counts/count_ind 

# mv counts/Matrix_data.tabular ~/

# ## Normalisation
# Rscript edgeR.r -m Matrix_data.tabular -a Galaxy_An/gffread_on_data_6__gtf.gtf -i Stress_level::NonStr,NonStr,NonStr,Str,Str,Str -o ./Galaxy_An -C Str,NonStr -z 1 -y
# Rscript deseq2.r -m Matrix_data.tabular -i NonStr,NonStr,NonStr,Str,Str,Str 
