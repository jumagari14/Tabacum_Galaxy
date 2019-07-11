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

# On terminal, paste sbatch -J RNA-Seq -N 16 --tasks-per-node 6 -t 20:00:00 ./main_rna_seq.sh

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
module load R/3.5.0 ## Load corresponding R module available
module load jdk1.8/8u22 ## Load corresponding java module available 

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib

## java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 /media/jumagari/JUANMA/Stage/Galaxy_An/subdata/TAB0.3_1.fq.gz /media/jumagari/JUANMA/Stage/Galaxy_An/subdata/TAB0.3_2.fq.gz forw_par.gq.gz forw_unp.fq.gz rev_pair.fq.gz rev_unp.fq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50
cd ~/ 
list=$(ls ./data/subdata/*.fq.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

mkdir -p -m 755 ./bin
mkdir -p -m 755 ./bin/trimm_data 
mkdir -p -m 755 ./bin/trimm_data/quality


gff=$(find ./data/ -not -path '*/\R*' -name "*gff[3]*")
fasta=$(find ./data/ -name "*.fa")

~/scripts/gffread "$gff" -o "$gff".gtf
tail -n +4 "$gff".gtf > ~/data/tmp ; mv ~/data/tmp "$gff".gtf

gtf=$(find ./data/ -not -path '*/\R*' -name "*.gtf")
gtf=$(readlink -f $gtf)
if [ ! -d "./bin/genome_ind" ] 
 then  ## Genome indexes are generated if necessary 
    mkdir -m 755 -p ./bin/genome_ind 
    ./scripts/STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ./bin/genome_ind  \
    --genomeFastaFiles "$fasta" \
    --sjdbOverhang 149 \
    --sjdbGTFfile "$gtf" \
    --genomeChrBinNbits 12 # reduce RAM consumption
fi  

## Only to be run once 
# if [$1=="convert"] 
# then 
#     parallel -j 3 "mv subdata/{}_1.fq.gz subdata/{}_1.fastqsanger" ::: $list 
#     parallel -j 3 "mv subdata/{}_2.fq.gz subdata/{}_2.fastqsanger" ::: $list
# fi

# Uncompress gz files 
parallel -j 2 "gunzip {}_1.fq.gz" ::: $list

mkdir -p -m 755 ./bin/STAR_Align  
mkdir -p -m 755 ./bin/counts

parallel -j 3 "./scripts/FastQC/fastqc ./data/subdata/{}_1.fq.gz ./data/subdata/{}_2.fq.gz --outdir=./bin/trimm_data/quality" ::: $list

parallel -j 3 "java -jar ./scripts/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 -phred33 ./data/subdata/{}_1.fq.gz ./data/subdata/{}_2.fq.gz ./bin/trimm_data/{}_1.par.fq.gz {}_1.unp.fq.gz ./bin/trimm_data/{}_2.par.fq.gz {}_2.unp.fq.gz ILLUMINACLIP:./scripts/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50" ::: $list

parallel -j 3 "rm -rf {}_1.unp.fq.gz {}_2.unp.fq.gz" ::: $list

./scripts/STAR --genomeLoad LoadAndExit --genomeDir ./bin/genome_ind
# Load genome just once to save RAM memory 

parallel --compress -j 3 "./scripts/STAR --runThreadN 16 --runMode alignReads --genomeLoad LoadAndKeep --readFilesCommand zcat --outFilterMatchNmin 16 --outSAMtype BAM Unsorted SortedByCoordinate --limitBAMsortRAM 10000200000 --genomeDir ./bin/genome_ind  --outFileNamePrefix ./bin/STAR_Align/{} --readFilesIn ./bin/trimm_data/{}_1.par.fq.gz  ./bin/trimm_data/{}_2.par.fq.gz" ::: $list
# --genomeLoad LoadAndKeep: keep genome in memory  
# --outFilterMatchNmin: alignment is set if matched bases is higher 
# --readFilesCommand zcat: work with compressed gz files
# --runThreadN: number of threads
# --runMode alignReads: alignment operation is run
# --outSAMtype BAM Unsorted SortedByCoordinate: both sorted and unsorted bam files are calculated
./scripts/STAR --genomeLoad Remove --genomeDir ./bin/genome_ind 
# Unload genome

# In featureCounts, paired-end reads must have -p option!!! 
cd ./bin/STAR_Align
parallel -j 3 "~/scripts/featureCounts -p -t exon -g gene_id -F GTF -Q 32 -T 16 -a "$gtf" -o ../counts/{.}.txt {}" ::: *Aligned.out.bam 
# p: mandatory option with paired-end reads
# -F GTF: reference format, GFF or GTF
# -t exon: exon region are kept from gtf file
# -Q 32: quality threshold
# -T 16: number of threads
# -a "$gtf": annotation file
# -o ../counts/{.}.txt: output prefix 


cd ../counts ; touch Matrix_data.tabular 

count=0

for J in $list 
do 
    count=$((count+1))
    if (($count == 1))
    then 
        sed '1d' "$J"Aligned.out.txt > tmp ; mv tmp "$J"Aligned.out.txt 
        cut -d $'\t' -f 1,7 "$J"Aligned.out.txt  > count_ind
        paste Matrix_data.tabular count_ind > temp2 && mv temp2 Matrix_data.tabular 
        rm -f temp2
    else
        sed '1d' "$J"Aligned.out.txt > tmp ; mv tmp "$J"Aligned.out.txt 
        cut -d $'\t' -f 7 "$J"Aligned.out.txt  > count_ind
        paste Matrix_data.tabular count_ind > temp2 && mv temp2 Matrix_data.tabular 
        rm -f temp2 
    fi 
done 

cd ../; rm -f counts/count_ind 

mv counts/Matrix_data.tabular ./


# ## Normalisation (to be run on the computer in a local folder)
# Rscript edgeR.r -m Matrix_data.tabular -a Galaxy_An/gffread_on_data_6__gtf.gtf -i Stress_level::Str,NonStr,Str,NonStr,Str,NonStr -o ./Galaxy_An -C Str,NonStr -z 1 -y -l 1
# -m: Path to count matrix
# -a: Path to input containing gene annotations
# -i: Factors with levels separated by comma 
# -o: Path to folder to write all output to
# -z: Integer specifying minimum total count requirement
# -y: Boolean to specify if total requirement is considered 
# -l: Float specifying the log-fold-change requirement   
# Rscript deseq2.r -m Matrix_data.tabular -i Str,NonStr,Str,NonStr,Str,NonStr -r 5
# -m: Path to count matrix
# -i: Levels separated by comma
# -r: Minimum count requirement per each feature