#/usr/local/bin/python3.9

import Bio.motifs as motifs
import pyfastx
import numpy as np
import pandas as pd
import os
print(os.chdir(os.path.dirname(__file__)))

# editing the text file to make it compatible for motif package.
data=pd.read_csv('argR-counts-matrix.txt',delimiter='\t',header=None)
S=data.iloc[0:4,2:20]# choose only numeric value for file.
array=S.to_numpy()# converted the selected data into numpy array.
np.savetxt("countmatrix.txt", array)# save the values in text file.

# opening the count matrix fie, to convert it into PSSM MATRIX
with open("countmatrix.txt") as handle:
 m = motifs.read(handle, "pfm")
pwm = m.counts.normalize(pseudocounts=1)# Added psudocounts and normalize the matrix.
background = {"A":0.25,"C":0.25,"G":0.25,"T":0.25}
pssm = pwm.log_odds(background)
print("The PSSM\n", pssm)
print("The consensus seq",pssm.consensus)#pssm matrix
print("The pssm max %4.2f" % pssm.max)#max value
print("The pssm min %4.2f" % pssm.min)#min value

#**************** to find out the threshold value to get relevant matches.
#******* it takes time to calculate threshold value so i just put the value after running this script.

#distribution = pssm.distribution(background=background, precision=10**4)
#threshold = distribution.threshold_balanced(1000)
#print("%5.3f" % threshold)
Threshold=8.508

#*************************
# Opening and bz2 file and coverted it into acceptable fasta format.
f=pd.read_csv('E_coli_K12_MG1655.400_50.bz2',delimiter=' ',header=None)
f.iloc[:,0] = '>' + f.iloc[:,0].astype(str)
f.drop(f.columns[1], axis=1, inplace=True)
f.drop(f.columns[2],axis=1,inplace=True)
f.to_csv('ecoli.fasta.txt',index=None,header=None)
with open(r'ecoli.fasta.txt', 'r') as file:
  data = file.read()
  data = data.replace(',', '\n')
with open(r'ecoli_formated.fasta', 'w') as file:
    file.write(data)
file.close()

# running the fasta file to get scores.
# saved the ouput in file scores.csv
fa=pyfastx.Fastx('ecoli_formated.fasta')
for name,seq,comment in fa:
    for position, score in pssm.search(seq, threshold=8.508):
        print("> %s:,Position %d: score = %5.3f" % (name,position, score),file=open("scores.csv", "a"))
file.close()

#opening the scores.csv file to get top 30 record ID for best transcription binding sites.

colnames=['RecordID', 'Position', 'scores',]
score_data= pd.read_csv('scores.csv',delimiter=':', names=colnames, header=None)
score_data['scores'] = score_data['scores'].map(lambda x: str(x)[8:])
score_data['Position'] = score_data['Position'].map(lambda x: str(x)[10:])
score_data['scores'] = score_data['scores'].apply(pd.to_numeric)
N=1
d=score_data.groupby('RecordID').head(N).reset_index(drop=True)
s=d.nlargest(30, 'scores')
s.to_csv("top30records.csv",index=False)
