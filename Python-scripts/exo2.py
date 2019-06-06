import argparse

from pathlib import Path

parser=argparse.ArgumentParser(prog="Gff file manipulator", description="This scripts extracts information from a gff file and stores it in a new one.")

parser.add_argument('-i',action='store',dest='input_file',type=Path,help='GFF filename that contains the data')
parser.add_argument('-o',action='store',dest='output_file',type=Path,help='Output GFF filename')
parser.add_argument('-m', action='store',dest='match',type=str,help='Match conditions')

args=parser.parse_args()



filename=args.output_file

f=open(args.input_file,'r')
new_f=open(filename,'w')

for x in f:
    line=x.strip('\n')
    line=line.split('\t')
    if (line[2]==args.match) : 
        last=line[len(line)-1].split(';') # A new list containing every feature from the last column separated by a semicolon is calculated
        search="Target="
        res=[i for i in last if search in i] # Information about 'Target' is extracted from the previously obtained list
        new_res=res[0].replace("Target=","") # Substring 'Target=' is removed
        for j in range(0,len(line)-1) :
            new_f.write(line[j]+"\t")
        new_f.write(new_res+"\n") # In the last column, information about 'Target' is just added

f.close()
new_f.close

