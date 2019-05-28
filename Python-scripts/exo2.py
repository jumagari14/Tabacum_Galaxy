import argparse

parser=argparse.ArgumentParser(prog="Gff file manipulator", description="This scripts extracts information from a gff file and stores it in a new one.")

parser.add_argument('-i',action='store',dest='input_file',type=str,help='GFF filename that contains the data')
parser.add_argument('-o',action='store',dest='output_file',type=str,help='Output GFF filename')
parser.add_argument('-m', action='store',dest='match',type=str,help='Match conditions')

args=parser.parse_args()



filename=args.output_file

f=open(args.input_file,'r')
new_f=open(filename,'w')

for x in f:
    line=x.strip('\n')
    line=line.split('\t')
    if (line[2]=='match') : 
        last=line[len(line)-1].split(';')
        search="Target="
        res=[i for i in last if search in i]
        new_res=res[0].replace("Target=","")
        for j in range(0,len(line)-1) :
            new_f.write(line[j]+"\t")
        new_f.write(new_res+"\n")

f.close()
new_f.close

