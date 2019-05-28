import argparse 


""" parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')

args = parser.parse_args()
print args.accumulate(args.integers) """

parser=argparse.ArgumentParser(prog="Little program", description="This is a little program")

parser.add_argument('-i',action='store',dest='input_file',type=str)
parser.add_argument('-o',action='store',dest='output_file',type=str)

args=parser.parse_args()

file = open(args.output_file,"w") 
 
file.write("Hello World\n") 
file.write("This is our new text file.\n") 
file.write("and this is another line.\n") 
file.write("Why? Because we can.") 
 
file.close() 