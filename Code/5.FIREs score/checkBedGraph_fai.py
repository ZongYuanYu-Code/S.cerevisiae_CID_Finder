import sys

bedgraph=open(sys.argv[1],'r')
fai=open(sys.argv[2],'r')
out=open(sys.argv[3],'w')

length={}
for line in fai:
	temp=line.rstrip().split()
	length[temp[0]]=int(temp[1])
	
for line in bedgraph:
	temp=line.rstrip().split()
	if int(temp[2])>length[temp[0]]:
		out.write(temp[0]+"\t"+temp[1]+"\t"+str(length[temp[0]])+"\t"+temp[3]+"\n")
	else:
		out.write(line)

