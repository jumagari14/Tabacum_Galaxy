Tabacum_Galaxy
==============
This repository stores all the Python, R and shell files that are used to perform NGS analysis.

## RNA-Seq folder
On the _RNA-Seq_ folder, you will find the necessary tools to run a RNA-Seq pipeline on your terminal in a Unix environment. Several folder are located within it: 
* data: _fasta_, _gff_ and a subfolder _subdata_ must be added in this folder. 

   * subdata: _fastq_ files are included in this folder. 
   
* scripts: all the scripts needed to run the _RNA-Seq_ analysis are found here. 

A folder _bin_ must be created where all the results will be stored. 
The user must follow the following steps before running the analysis 
ºººsh 
git clone https://github.com/jumagari14/Tabacum_Galaxy.git
cd Tabacum_Galaxy/RNA-Seq
mkdir -p -m 755 bin
cd scripts 
./improved-rna.alignment.sh 
