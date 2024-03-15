import sys
import os
import argparse
import fire

def modify(fai, border, res):
	res = int(res)
	with open(fai, 'r')as f:
		for line in f:
			line = line.rstrip().split()
			chr_id = line[0]
			chr_len = int(line[1])
			i = 0
			z = open(str(chr_id)+"-"+str(res)+".border.bed", 'w')
			with open(border, 'r')as g:
				for info in g:
					info = info.rstrip().split()
					if info[0] == chr_id:
						i += 1
						border_bin = int(info[1])/res + 1
						z.write(str(chr_id)+"_"+str(i)+'\t'+str(border_bin)+'\n')
	
			z.close()


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fai", type=str, help="the chr length file")
	parser.add_argument("-r", "--res", type=int, help="the resolution of the fire file")
	parser.add_argument("-d", "--border", type=str, help="the tad border file")
	command = sys.argv[1:]
	args = parser.parse_args(command)
	return command,args

def main():
	command,args = get_args()
	fai = args.fai
	res = args.res
	border = args.border
	modify(fai, border, res)	

if __name__ == "__main__":
	 main()
