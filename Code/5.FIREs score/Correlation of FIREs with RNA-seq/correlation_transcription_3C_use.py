# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:31:33 2017
@author: Axel KournaK
To look for correlation between transcription and 3C contacts. 
"""
import sys
#sys.path.insert(0,'/public/frasergen/3D/pipeline/Interactome/bacteria/E_coli_analysis_tools_2018_Cell/python_codes')
#from pylab import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
#import ice_mirny3
import scn
from itertools import islice
from scipy.stats import gaussian_kde
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from scalogram_2d_type import *
import sys
# smooth function:
def smoothListGaussian(list,strippedXs=False,degree=5):  
    window=degree*2-1  
    weight=np.array([1.0]*window)  
    weightGauss=[]  
    for i in range(window):  
        i=i-degree+1  
        frac=i/float(window)  
        gauss=1/(np.exp((4*(frac))**2))  
        weightGauss.append(gauss)  
    weight=np.array(weightGauss)*weight  
    smoothed=[0.0]*(len(list)-window)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  
    return smoothed


data1=np.loadtxt(sys.argv[1])

#matscn1=ice_mirny3.ice_func(data1,100)
#matscn1= scn.scn_func(matscn1,0)
b=data1[:,1]
print (b)
BIN= 1000
# plots:

#b=np.diagonal(matscn1)
#b=scalo2(matscn1,[2],BIN);
#print(b)
#plot(b)
#plot(smoothListGaussian(b) )

# transcription data (Rnaseq from Olivier Espeli lab):
#chip=loadtxt("/run/media/axel/9e657c5d-6ac3-494e-81af-b25e389d59bd1/vick_data_backup/espeli_data/fastq/EV-4_TGACCA_L002_R1_001.fastq.sam.MQ0.hist5000")
chip=np.loadtxt(sys.argv[2])
#plot(log(chip))

# Correlation:
c=np.log(chip[:,1])
print (c)
#plot( (b-mean(b))/std(b),color="Blue")
#plot( (c-mean(c))/std(c),color="Green")

ax=pearsonr(b,c)
num=len(c)
print("pearsonr")
print(ax)
print (spearmanr(b, c))
# with smoothing the signals:
bx=pearsonr(smoothListGaussian( (b-np.mean(b))/np.std(b)  ), smoothListGaussian( (c-np.mean(c)/np.std(c) )  ))
print ("aaa")
cx=spearmanr(smoothListGaussian( (b-np.mean(b))/np.std(b)  ), smoothListGaussian( (c-np.mean(c)/np.std(c) )  ))
print ("bbb")
#bx=pearsonr((b-mean(b))/std(b), (c-mean(c)/std(c)))
#cx=spearmanr((b-mean(b))/std(b),(c-mean(c)/std(c)))
print("pearsonr_smoothListGaussian")
print(bx)
print("spearmanr_smoothListGaussian")
print(cx)
print("#" * 100)
print(os.environ.get('DISPLAY'))
plt.figure(figsize=(8, 2))
#  plot of both signal:
x_zhou=range(5,len(smoothListGaussian((b-np.mean(b))/np.std(b)))+5)
plt.plot(x_zhou,smoothListGaussian( (b-np.mean(b))/np.std(b) ),color="DodgerBlue",linewidth=2,label = "FIRE signal")
plt.plot(x_zhou,smoothListGaussian( (c-np.mean(c))/np.std(c) ),color="Red",linewidth=2,label = "transcription binned at {}kb".format(sys.argv[5]))
#plt.plot( (b-mean(b))/std(b),color="DodgerBlue",linewidth=2,label = "FIRE signal")
#plt.plot((c-mean(c))/std(c),color="Red",linewidth=2,label = "transcription binned at 5kb")
plt.title("PC:{:.4f},p-value:{:2e}".format(bx[0],bx[1]))
ylabel_obj=plt.ylabel('Z transformed signal',fontsize=10)
plt.xlim(0, num)
#plt.text(100,2.5,"PC:{:.4f}\np-value:{:2e}".format(bx[0],bx[1]))
with open(sys.argv[3],'r') as border:
  for line in border:
    if (line[0:11] == 'boundary_id'):
      continue
    line_split=line.strip().split()
    border_bin=int(float(line_split[1]))
    linemid=float(border_bin)-0.5
    plt.axvline(linemid,color="Gray",linewidth=0.5,linestyle='--')
plt.legend(loc = 0, ncol = 1)
#show()
out=sys.argv[4]
plt.savefig(out);

