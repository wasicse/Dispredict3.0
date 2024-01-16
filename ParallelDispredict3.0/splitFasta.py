from Bio import SeqIO
import sys
import math
import pandas as pd
print ('Number of parallel runs: ', sys.argv[1])
# Number of sequences in each file

i=1
count1=0

df=pd.DataFrame(columns=['Protein','Sequence'])
with open('../example/sample.fasta', 'r') as f:
  for record in SeqIO.parse(f, 'fasta'):
        df2=pd.DataFrame({'Protein':[str(record.id)],'Sequence':[str(record.seq)]})
        df=pd.concat([df,df2],ignore_index=True)
        count1 += 1
        
df = df.sample(frac = 1)
        
print ('Total number proteins: ', count1)
splitcount=math.ceil(count1/int(sys.argv[1]))
print('Number of proteins in each file: ',splitcount)


print(df.shape)

count=0  
i=1

for p in range(count1):

        
        if(count<=splitcount):
          # print(record.id)
          count += 1
          newfile=open('./temp/Parallelinputs/processedinput_'+str(i)+'.fasta','a')
    
          newfile.write('>' + df.iloc[p,0] + '\n')
          newfile.write(str(df.iloc[p,1] ) + '\n')
        if(count==int(splitcount)):
          i=i+1
          
          count=0
          newfile.close()
if not newfile.closed: 
  newfile.close()    
i=i-1
print("Number of files", i)