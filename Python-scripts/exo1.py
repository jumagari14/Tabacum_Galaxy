""" import pprint
from BCBio.GFF import GFFExaminer

ex=GFFExaminer()
reader=open('reswag5S_REPET_SSRs.gff')

pprint.pprint(ex.available_limits(reader))

reader.close() """

import argparse

from pathlib import Path

parser=argparse.ArgumentParser(prog="Gff file manipulator", description="This scripts extracts information from a gff file and stores it in a new one.")

parser.add_argument('-i',action='store',dest='input_file',type=Path,help='GFF filename that contains the data')
parser.add_argument('-o',action='store',dest='output_file',type=Path,help='Output GFF filename')
parser.add_argument('-m', action='store',dest='match',type=str,help='Match conditions')
parser.add_argument('-t', action='store',nargs='+',dest='target',type=str,help='Target conditions, separated by a space. If space separator are included in the conditions, the escape character \ must be used')

args=parser.parse_args()

filename=args.output_file

f=open(args.input_file,'r')
new_f=open(filename,'w')


for x in f:
    line=x.strip('\n')
    line=line.split('\t') # In one row, a list containing data from one column in each position is generated
    for i in range(0,len(args.target)) :         
        if (line[2]==args.match and ("Target=" + args.target[i] in line[len(line)-1])) : # In order to write on the new file, the substring 'Target=input target from terminal ' is searched in the last element of the list, which corresponds to the last column
                new_f.write(x)
f.close()
new_f.close

        
   

