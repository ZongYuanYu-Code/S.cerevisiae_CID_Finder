from optparse import OptionParser
#from mirnylib.h5dict import h5dict
import numpy
import os
import sys
def parseArgs():
    parser = OptionParser()
    parser.add_option("-b", "--binRange", dest="binRange", type="int", default=200000,\
                        help="Bin Range for creating local interaction file. ")
    parser.add_option("-i", "--input", dest="input",\
                        help="Input ICE generated raw hm file")
    parser.add_option("-r", "--res", dest="res", type="int", default=0,\
                        help="resolution of the heatmap. DEFAULT=0")
    parser.add_option("-c", "--characterFile",dest="characterFile", type=str,\
                        help=("the characterFile file of the genome" ))

    parser.add_option("-o", "--out",dest="out", type=str,\
                        help=("the characterFile file of the genome" ))
    parser.add_option("-a", "--fai",dest="fai", type=str,\
                        help=("the fai file of the genome" ))

    #Put this in. I think it would be nice.
    #parser.add_option("-o","--output", dest="output", type="str", default=None,
                        #help="File name to output a Bed File with chromosome, start, end of gaps features. Default=Off")
    opts,args = parser.parse_args()
    '''
    if len(args) != 1:
        parser.error("No Fasta Specified!")
    if opts.min < 1:op
        parser.error("Minimum gap size must be at least 1.")
    '''
    return opts,args	
	


#def get_chrInfo(hm_file):
#	cmd="/public/frasergen/PUB/software/Anaconda/anaconda3-3d/envs/iced/bin/python /public/frasergen/3D/pipeline/Interactome/animal/ice_pipe/src/Fire/h5dictToTxt.py "+hm_file + " temp"
#	os.system(cmd)
#	return os.path.abspath("temp/chromosomeIndex.txt")

def hm2FIRE(hm_file,res,binNum,chrFile): ## build the local interaciton 3 colum file
	#dataset = h5dict(hm_file,'r')
	#print 'About to loading heatmap!!'
	#chromosomeHeatmap = dataset["heatmap"]
	print 'About to loading heatmap!!'
	chromosomeHeatmap=numpy.loadtxt(hm_file)
	print 'loading completed!!!'
	arraySize=chromosomeHeatmap.shape[0]

#chrInfo=open(chrFile,'r')
	chrN=numpy.loadtxt(chrFile)
	#outfilename=os.path.basename(hm_file).split(".hm")[0]+"_ForFIRE.txt"
	outfilename=os.path.basename(hm_file).split(".matrix")[0]+"_ForFIRE.txt"
	out=open(outfilename,'w')
	i=0
	QianBin=0
	QianChr=''
	dangqian=0
	#print binNum
	while i<arraySize:
		if i<int(binNum):
			tt=numpy.sum(chromosomeHeatmap[0:i+int(binNum)+1,i])
			oo="chr"+str(int(chrN[i]+1))+"\t"+str(i*int(res))+"\t"+str(tt)+"\n"
			QianChr="chr"+str(int(chrN[i]+1))
			out.write(oo)
		else:
			tt=numpy.sum(chromosomeHeatmap[i-int(binNum):i+int(binNum)+1,i])
		#print tt
			if "chr"+str(int(chrN[i]+1))==QianChr:
				a=dangqian*int(res)
				oo="chr"+str(int(chrN[i]+1))+"\t"+str(a)+"\t"+str(tt)+"\n"
				out.write(oo)
			else:
				dangqian=0
				a=dangqian*int(res)
				oo="chr"+str(int(chrN[i]+1))+"\t"+str(a)+"\t"+str(tt)+"\n"
				out.write(oo)
		#print 'aa'
		QianChr="chr"+str(int(chrN[i]+1))
		i+=1
		dangqian+=1
	return outfilename

def formatOut(out,fai,res):
	cmd1="awk '{print $1,$2,$2+"+str(res)+",$3} ' "+ out+" |grep -v 'NaN'> "+out+".bedgraph"
	cmd2="python checkBedGraph_fai.py "+out+".bedgraph "  +fai+ " "+out+".bedgraph.checked"
	cmd3="awk '$4>2' "+out+".bedgraph.checked >"+out+".bedgraph.checked.FIRE"
	print cmd1
	print cmd2
	print cmd3
	os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	
if __name__ == '__main__':
	opts,args = parseArgs()
#	print opts["input"]
	inputHmFile=os.path.abspath(opts.input)
	#chrFile=get_chrInfo(inputHmFile)
	chrFile="/mnt/g/S.cere-work/8.14-HiCNormCis/chromosomeIndex5kb.txt"
	binNum=opts.binRange/opts.res
	localInterFile=hm2FIRE(inputHmFile,opts.res,binNum,chrFile)
	cmd1="Rscript /mnt/g/S.cere-work/8.14-HiCNormCis/HiCNormCis.R -i "+localInterFile+" -f "+opts.characterFile +" -o   "+opts.out
	os.system("Rscript /mnt/g/S.cere-work/8.14-HiCNormCis/HiCNormCis.R -i "+localInterFile+" -f "+opts.characterFile +" -o   "+opts.out)
	formatOut(opts.out,opts.fai,opts.res)	


