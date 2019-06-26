Tabacum_Galaxy
==============
This repository stores all the Python, R and shell files that are used to perform NGS analysis.

## RNA-Seq pipeline
On the _RNA-Seq_ folder, you will find the necessary tools to run a RNA-Seq pipeline on your terminal in a Unix environment. Several folder are located within it: 
* data: _fasta_, _gff_ and a subfolder _subdata_ must be added in this folder. 

   * subdata: _fastq_ files are included in this folder. 
   
* scripts: all the scripts needed to run the _RNA-Seq_ analysis are found here. 

A folder _bin_ must be created where all the results will be stored. 
The user must follow the following steps before running the analysis 
```bash 
git clone https://github.com/jumagari14/Tabacum_Galaxy.git
cd Tabacum_Galaxy
tar xvzf RNA-Seq.tar.gz 
cd RNA-Seq
mkdir -p -m 755 bin data
mkdir -p -m 755 data/subdata
## Include all the necessary data in data and subdata folder
cd scripts 
```
Once this steps are done, a _shell_ file can be executed 
```bash
./improved_rna_seq.sh 
```

After running the pipeline, several folders are found in _bin_: 
* genome_ind: Genome indexes, necessary to run the STAR mapping. 
* STAR_Align: Sorted and unsorted _Bam_ files from STAR mapping.
* counts: _txt_ files with the counting results. 
* trimm_data: trimmed _fastq_ reads. A subfolder _quality_ where _html_ and _zip_ files as a result of _FASTQC_ analysis are stored is also generated. 

In the main directory, a file _.tabular_ where all the counting results are saved is created. This file will mainly be used in the latter normalisation step.  
