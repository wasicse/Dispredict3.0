#!/home/pyenv/versions/py3.7/bin/python
# coding: utf-8
import numpy as np
import pandas as pd
import io
from io import StringIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import os
# In[2]:

warning_code = "[WARNING]"
warning_msg =  " [WARNING-DEFAULT PSSM]"

input_csv_path  = sys.argv[1]
output_dir = os.path.dirname(input_csv_path)



CAID_ID=[]
Length_EachProtein = []
with open(input_csv_path) as resf:
	spl = resf.read().split(">")
	lines = spl[1:][0].split("\n")
	seqid = lines[0]
	CAID_ID.append(seqid)
	df_str = "\n".join(lines[1:])
	Prediction_Df = pd.read_csv(StringIO(df_str))
	Length_EachProtein.append(len(Prediction_Df))


# In[3]:


b=1
for b in range(1,len(spl[1:]),1):
	lines = spl[1:][b].split("\n")
	seqid = lines[0]
	CAID_ID.append(seqid)
	df_str = "\n".join(lines[1:])
	Df = pd.read_csv(StringIO(df_str))
	Length_EachProtein.append(len(Df))
	Prediction_Df =[Prediction_Df,Df]
	Prediction_Df= pd.concat(Prediction_Df)


# In[4]:


d=np.cumsum(Length_EachProtein)


# In[5]:


Prediction_Df.head()


# In[6]:


rdp_r = Prediction_Df["rdp_r"].values
rdp_d = Prediction_Df["rdp_d"].values
rdp_p = Prediction_Df["rdp_p"].values
dfl = Prediction_Df["dfl"].values
fmorf = Prediction_Df["fmorf"].values
AA_Sequences= Prediction_Df["Residue Type"].values
DisorderPrediction_Score=Prediction_Df["Predicted Score for Disorder"].values


# In[7]:


DisorderBinary_Prediction=[]
b=0
for b in range(0,len(DisorderPrediction_Score),1):
	if DisorderPrediction_Score[b]>0.3:DisorderBinary_Prediction.append(1)
	else:DisorderBinary_Prediction.append(0)


# In[8]:


DisorderPrediction_Scores= np.split(DisorderBinary_Prediction,d)
ModifiedBinary_Predictions= np.split(DisorderBinary_Prediction,d)


rdp_r_byProteins= np.split(rdp_r,d)
rdp_d_byProteins= np.split(rdp_d,d)
rdp_p_byProteins= np.split(rdp_p,d)
dfl_byProteins= np.split(dfl,d)
fmorf_byProteins= np.split(fmorf,d)

AA_Sequences_byProteins= np.split(AA_Sequences,d)


# In[9]:


def window_calc(NewScore10_byProteins):
	Prediction_ScoreWindow= []
	newscore10window=0
	b = 0
	for NewScore10_byProtein in NewScore10_byProteins:
		for b in range(0, len(NewScore10_byProtein), 1):
			if b == 0 : newscore10window= (NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/6
			elif b == 1 : newscore10window= (NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/7
			elif b == 2 : newscore10window= (NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/8
			elif b == 3 : newscore10window= (NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/9
			elif b == 4 : newscore10window= (NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/10
			elif len(NewScore10_byProtein)-b == 5 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4])/10
			elif len(NewScore10_byProtein)-b == 4 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3])/9
			elif len(NewScore10_byProtein)-b == 3 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2])/8
			elif len(NewScore10_byProtein)-b == 2 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1])/7
			elif len(NewScore10_byProtein)-b == 1 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b])/6
			else: newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/11
			Prediction_ScoreWindow.append(newscore10window)

	Prediction_ScoreWindowWeighted= []
	newscore10windowweighted=0
	b = 0
	for NewScore10_byProtein in NewScore10_byProteins:
		for b in range(0, len(NewScore10_byProtein), 1):
			if b == 0 : newscore10windowweighted= (NewScore10_byProtein[b]*0.65+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/6
			elif b == 1 : newscore10windowweighted= (NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.56+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/7
			elif b == 2 : newscore10windowweighted= (NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.48+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/8
			elif b == 3 : newscore10windowweighted= (NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.41+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/9
			elif b == 4 : newscore10windowweighted= (NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.35+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/10
			elif len(NewScore10_byProtein)-b == 5 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.35+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06)/10
			elif len(NewScore10_byProtein)-b == 4 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.41+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07)/9
			elif len(NewScore10_byProtein)-b == 3 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.48+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08)/8
			elif len(NewScore10_byProtein)-b == 2 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.56+NewScore10_byProtein[b+1]*0.09)/7
			elif len(NewScore10_byProtein)-b == 1 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.65)/6
			else: newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.3+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/11
			Prediction_ScoreWindowWeighted.append(newscore10windowweighted)
	return Prediction_ScoreWindow, Prediction_ScoreWindowWeighted


# In[10]:


rdp_r_ScoreWindow,rdp_r_ScoreWindowWeighted=window_calc(rdp_r_byProteins)
rdp_d_ScoreWindow,rdp_d_ScoreWindowWeighted=window_calc(rdp_d_byProteins)
rdp_p_ScoreWindow,rdp_p_ScoreWindowWeighted=window_calc(rdp_p_byProteins)
dfl_ScoreWindow,dfl_ScoreWindowWeighted=window_calc(dfl_byProteins)
fmorf_ScoreWindow,fmorf_ScoreWindowWeighted=window_calc(fmorf_byProteins)
Disorder_ScoreWindow,Disorder_ScoreWindowWeighted=window_calc(DisorderPrediction_Scores)


# In[11]:


#######################################################################################################
#Getting the positions of predicted Disorder regions
Disorder_Index= []
disorderindex=0
b = 0
for b in range(0, len(DisorderBinary_Prediction), 1):
	if DisorderBinary_Prediction[b] == 1 : disorderindex= b
	Disorder_Index.append(disorderindex)
Disorder_Index=np.unique(Disorder_Index)

######################################################################################################
#Getting the margins of disorder regions
Disorder_Margins= []
disordermargins=0
b = 0
h = len(Disorder_Index) -1
for b in range(0, h , 1):
	if Disorder_Index[b+1]-Disorder_Index[b] != 1 : disordermargins= (b+1)
	elif np.isin(Disorder_Index[b], d) == True : disordermargins= (b)
	Disorder_Margins.append(disordermargins)
Disorder_Margins=np.unique(Disorder_Margins)

##############################################################################################################
#Splitting the disorder positions in to regions
Disorder_Regions= np.split(Disorder_Index,Disorder_Margins)

##############################################################################################################
#Calculating the depth of all disorder residues in their regons
Disorder_Depth= []
disorderdepth=0
b = 0
for Disorder_Region in Disorder_Regions:
	for b in range(0, len(Disorder_Region), 1):
		if b < (len(Disorder_Region))/2 : disorderdepth = b
		elif b > (len(Disorder_Region))/2 : disorderdepth = len(Disorder_Region) -(b+1)
		Disorder_Depth.append(disorderdepth)

################################################################################################################
#Assigning -1 to ordered residues
Disorder_DepthScore= []
disorderdepthscore=0
b = 0
for b in range(0, len(DisorderBinary_Prediction), 1):
	if DisorderBinary_Prediction[b] == 0 : disorderdepthscore= -1
	else: disorderdepthscore= DisorderBinary_Prediction[b]
	Disorder_DepthScore.append(disorderdepthscore)

#############################################################################################################
#Casting lists to arrays
Disorder_DepthScore = np.asarray(Disorder_DepthScore)
Disorder_Index = np.asarray(Disorder_Index)
Disorder_Depth = np.asarray(Disorder_Depth)

#Inserting depths to the position of disorder residues
np.put(Disorder_DepthScore, Disorder_Index, Disorder_Depth)
Disorder_DepthScore=list(Disorder_DepthScore)


# In[12]:


#######################################################################################################
#Getting the positions of predicted Order regions
Order_Index= []
orderindex=0
b = 0
for b in range(0, len(DisorderBinary_Prediction), 1):
	if DisorderBinary_Prediction[b] == 0 : orderindex= b
	Order_Index.append(orderindex)
Order_Index=np.unique(Order_Index)

######################################################################################################
#Getting the margins of disorder regions
Order_Margins= []
ordermargins=0
b = 0
h = len(Order_Index) -1
for b in range(0,h , 1):
	if Order_Index[b+1]-Order_Index[b] != 1 : ordermargins= (b+1)
	elif np.isin(Order_Index[b], d) == True : ordermargins= (b)
	Order_Margins.append(ordermargins)
Order_Margins=np.unique(Order_Margins)

#############################################################################################################
#Splitting the disorder positions in to regions
Order_Regions= np.split(Order_Index,Order_Margins)

#############################################################################################################
#Calculating the depth of all disorder residues in their regons
Order_Depth= []
orderdepth=0
b = 0
for Order_Region in Order_Regions:
	for b in range(0, len(Order_Region), 1):
		if b < (len(Order_Region))/2 : orderdepth = b
		elif b > (len(Order_Region))/2 : orderdepth = len(Order_Region) -(b+1)
		Order_Depth.append(orderdepth)

################################################################################################################
#Assigning -1 to disordered residues
Order_DepthScore= []
orderdepthscore=0
b = 0
for b in range(0, len(DisorderBinary_Prediction), 1):
		if DisorderBinary_Prediction[b] == 1 : orderdepthscore= -1
		else: orderdepthscore= DisorderBinary_Prediction[b]
		Order_DepthScore.append(orderdepthscore)

#############################################################################################################
#Casting lists to arrays
Order_DepthScore = np.asarray(Order_DepthScore)
Order_Index = np.asarray(Order_Index)
Order_Depth = np.asarray(Order_Depth)

#Inserting depths to the position of disorder residues
np.put(Order_DepthScore, Order_Index, Order_Depth)
Order_DepthScore=list(Order_DepthScore)


# In[13]:


Terminal_Distance= []
terminaldistance=0
b = 0
for Prediction_Score_byProtein in DisorderPrediction_Scores:
	for b in range(0, (len(Prediction_Score_byProtein)), 1):
		if b <len(Prediction_Score_byProtein)/2 : terminaldistance= b/len(Prediction_Score_byProtein)
		else : terminaldistance = ((len(Prediction_Score_byProtein))-(b+1))/len(Prediction_Score_byProtein)
		Terminal_Distance.append(terminaldistance)


# In[14]:


Terminal_Distance10= []
terminaldistance10=0
b = 0
for Prediction_Score_byProtein in DisorderPrediction_Scores:
	for b in range(0, (len(Prediction_Score_byProtein)), 1):
		if b <10 : terminaldistance10= b
		elif (len(Prediction_Score_byProtein))-(b+1)<10 :terminaldistance10 = ((len(Prediction_Score_byProtein))-(b+1))
		else : terminaldistance10 = 10
		Terminal_Distance10.append(terminaldistance10)


# In[15]:


Feature_Set=pd.DataFrame({
'Binary_Prediction':DisorderBinary_Prediction,

'Disorder_Score':DisorderPrediction_Score,
'Disorder_ScoreWindow':Disorder_ScoreWindow,
'Disorder_ScoreWindowWeighted':Disorder_ScoreWindowWeighted,

'rdp_r_Score':rdp_r,
'rdp_r_ScoreWindow':rdp_r_ScoreWindow,
'rdp_r_ScoreWindowWeighted':rdp_r_ScoreWindowWeighted,

'rdp_d_Score':rdp_d,
'rdp_d_ScoreWindow':rdp_d_ScoreWindow,
'rdp_d_ScoreWindowWeighted':rdp_d_ScoreWindowWeighted,

'rdp_p_Score':rdp_p,
'rdp_p_ScoreWindow':rdp_p_ScoreWindow,
'rdp_p_ScoreWindowWeighted':rdp_p_ScoreWindowWeighted,

'dfl_Score':dfl,
'dfl_ScoreWindow':dfl_ScoreWindow,
'dfl_ScoreWindowWeighted':dfl_ScoreWindowWeighted,

'fmorf_Score':fmorf,
'fmorf_ScoreWindow':fmorf_ScoreWindow,
'fmorf_ScoreWindowWeighted':fmorf_ScoreWindowWeighted,

'Disorder_DepthScore':Disorder_DepthScore,
'Order_DepthScore':Order_DepthScore,
'Terminal_Distance':Terminal_Distance,
'Terminal_Distance10':Terminal_Distance10})
Feature_Set.head()


# In[16]:


Test_disorderindexes=[]
Test_orderindexes=[]
Test_orderscores=[]
b=0
for b in range(0,len(DisorderBinary_Prediction),1):
	if DisorderBinary_Prediction[b]==1:
		Test_disorderindexes.append(b)
	else:
		Test_orderindexes.append(b)
		Test_orderscores.append(DisorderPrediction_Score[b])


# In[17]:


Feature_Set=Feature_Set.iloc[Test_disorderindexes,:]


# In[18]:


X_test= np.array(Feature_Set.iloc[:,0:23])


# In[19]:


import pickle
clfProtein = pickle.load(open('Protein_Residue_DisoOnly_CV.sav','rb'))
clfDNA = pickle.load(open('DNA_Residue_DisoOnly_CV.sav','rb'))
clfRNA = pickle.load(open('RNA_Residue_DisoOnly_CV.sav','rb'))
clfLinker = pickle.load(open('Linker_Residue_DisoOnly_CV.sav','rb'))


Test_FullProteinPrediction = np.zeros((len(DisorderBinary_Prediction),),dtype=float)
Test_FullDNAPrediction = np.zeros((len(DisorderBinary_Prediction),),dtype=float)
Test_FullRNAPrediction = np.zeros((len(DisorderBinary_Prediction),),dtype=float)
Test_FullLinkerPrediction = np.zeros((len(DisorderBinary_Prediction),),dtype=float)

try:
	Test_ProteinScore = clfProtein.predict_proba(X_test)
	Test_DNAScore = clfDNA.predict_proba(X_test)
	Test_RNAScore = clfRNA.predict_proba(X_test)
	Test_LinkerScore = clfLinker.predict_proba(X_test)



	Test_ProteinScore=Test_ProteinScore[:,1].reshape(Test_ProteinScore[:,1].shape[0],1)
	Test_DNAScore=Test_DNAScore[:,1].reshape(Test_DNAScore[:,1].shape[0],1)
	Test_RNAScore=Test_RNAScore[:,1].reshape(Test_RNAScore[:,1].shape[0],1)
	Test_LinkerScore=Test_LinkerScore[:,1].reshape(Test_LinkerScore[:,1].shape[0],1)



	from sklearn.preprocessing import MinMaxScaler
	scaler = MinMaxScaler(feature_range=(0.3,1.0))
	scaler.fit(Test_ProteinScore)
	Test_ProteinScore=scaler.transform(Test_ProteinScore)
	Test_ProteinScore=np.concatenate(Test_ProteinScore)
	np.put(Test_FullProteinPrediction ,Test_disorderindexes,Test_ProteinScore)


	scaler = MinMaxScaler(feature_range=(0.3,1.0))
	scaler.fit(Test_DNAScore)
	Test_DNAScore=scaler.transform(Test_DNAScore)
	Test_DNAScore=np.concatenate(Test_DNAScore)
	np.put(Test_FullDNAPrediction ,Test_disorderindexes,Test_DNAScore)


	scaler = MinMaxScaler(feature_range=(0.3,1.0))
	scaler.fit(Test_RNAScore)
	Test_RNAScore=scaler.transform(Test_RNAScore)
	Test_RNAScore=np.concatenate(Test_RNAScore)
	np.put(Test_FullRNAPrediction ,Test_disorderindexes,Test_RNAScore)
	

	scaler = MinMaxScaler(feature_range=(0.3,1.0))
	scaler.fit(Test_LinkerScore)
	Test_LinkerScore=scaler.transform(Test_LinkerScore)
	Test_LinkerScore=np.concatenate(Test_LinkerScore)
	np.put(Test_FullLinkerPrediction ,Test_disorderindexes,Test_LinkerScore)

except:
	k=1

# In[27]:




np.put(Test_FullProteinPrediction ,Test_orderindexes,Test_orderscores)


# In[28]:



np.put(Test_FullDNAPrediction ,Test_orderindexes,Test_orderscores)


# In[29]:



np.put(Test_FullRNAPrediction ,Test_orderindexes,Test_orderscores)


# In[30]:



np.put(Test_FullLinkerPrediction ,Test_orderindexes,Test_orderscores)


# In[31]:




DisorderPrediction_Score=np.around(DisorderPrediction_Score,decimals=3)
Test_FullProteinPrediction=np.around(Test_FullProteinPrediction,decimals=3)
Test_FullDNAPrediction=np.around(Test_FullDNAPrediction,decimals=3)
Test_FullRNAPrediction=np.around(Test_FullRNAPrediction,decimals=3)
Test_FullLinkerPrediction=np.around(Test_FullLinkerPrediction,decimals=3)

Protein_Binary=[]
b=0
for b in range(0,len(Test_FullProteinPrediction),1):
	if Test_FullProteinPrediction[b]>0.5:Protein_Binary.append(1)
	else:Protein_Binary.append(0)


# In[32]:


DNA_Binary=[]
b=0
for b in range(0,len(Test_FullDNAPrediction),1):
	if Test_FullDNAPrediction[b]>0.5:DNA_Binary.append(1)
	else:DNA_Binary.append(0)


# In[33]:


RNA_Binary=[]
b=0
for b in range(0,len(Test_FullRNAPrediction),1):
	if Test_FullRNAPrediction[b]>0.5:RNA_Binary.append(1)
	else:RNA_Binary.append(0)


# In[34]:


Linker_Binary=[]
b=0
for b in range(0,len(Test_FullLinkerPrediction),1):
	if Test_FullLinkerPrediction[b]>0.5:Linker_Binary.append(1)
	else:Linker_Binary.append(0)


# In[35]:


DisorderPrediction_Score_byProteins= np.split(DisorderPrediction_Score,d)

ProteinPrediction_Score_byProteins= np.split(Test_FullProteinPrediction,d)
DNAPrediction_Score_byProteins= np.split(Test_FullDNAPrediction,d)
RNAPrediction_Score_byProteins= np.split(Test_FullRNAPrediction,d)
LinkerPrediction_Score_byProteins= np.split(Test_FullLinkerPrediction,d)


# In[36]:


DisorderBinary_Prediction_byProteins= np.split(DisorderBinary_Prediction,d)

ProteinBinary_Prediction_byProteins= np.split(Protein_Binary,d)
DNABinary_Prediction_byProteins= np.split(DNA_Binary,d)
RNABinary_Prediction_byProteins= np.split(RNA_Binary,d)
LinkerBinary_Prediction_byProteins= np.split(Linker_Binary,d)


# In[38]:

marker_size=8
j=0
#len(CAID_ID)
for j in range(0,len(CAID_ID),1):
	try:
			Disordered_Residues= np.concatenate(np.argwhere(DisorderBinary_Prediction_byProteins[j]>0))
	except:
			Disordered_Residues= list(np.argwhere(DisorderBinary_Prediction_byProteins[j]>0))
	Position =np.arange(0,len(DisorderBinary_Prediction_byProteins[j]),1)
	AminoAcids=list(AA_Sequences_byProteins[j])
	Protein_Score=ProteinPrediction_Score_byProteins[j]
	DNA_Score=DNAPrediction_Score_byProteins[j]
	RNA_Score=RNAPrediction_Score_byProteins[j]
	Linker_Score=LinkerPrediction_Score_byProteins[j]

	Protein_DisordredScores=[]
	DNA_DisordredScores=[]
	RNA_DisordredScores=[]
	Linker_DisordredScores=[]
	b=0
	for b in Disordered_Residues:
		Protein_DisordredScores.append(Protein_Score[b])
		DNA_DisordredScores.append(DNA_Score[b])
		RNA_DisordredScores.append(RNA_Score[b])
		Linker_DisordredScores.append(Linker_Score[b])

	Y_Sequence = np.full(len(Position), 6)
	Y_Disorder = np.full(len(Position), 5)
	Y_Protein = np.full(len(Position), 4)
	Y_DNA = np.full(len(Position), 3.5)
	Y_RNA = np.full(len(Position), 3)
	Y_Linker = np.full(len(Position), 2.5)

	buffer = io.StringIO()

	fig = make_subplots(rows=2, cols=1,
											shared_xaxes=True,
											vertical_spacing=0.01)



	fig.add_trace(go.Scatter(x=Position+1, y=DisorderPrediction_Score_byProteins[j],name='Disorder Score',
					marker=dict(color='black'),
					meta=AminoAcids,
					hovertemplate='Score: %{y:.3f}'+'<br>Position: %{x}'+'%{meta}'),row=1, col=1)

	fig.add_shape(
					type='line',
					x0=1,
					y0=0.3,
					x1=len(Position),
					y1=0.3,
					line=dict(
							color='black',
							dash="dot",
							width=3),
								row=1, col=1)

	fig.add_trace(go.Scatter(
			x=Position+1,
			y=Y_Sequence,
			mode="text",
			name="Sequence",
			text=AminoAcids,
			textposition="bottom center",
			marker=dict(size=2,color='black'),
			meta=AminoAcids,
			hovertemplate='Residue: %{meta}'+'<br>Position: %{x}'),
			 row=2, col=1 )

	fig.add_trace(go.Scatter(x=Disordered_Residues+1, y=Y_Disorder,mode='markers',name='Disordered Regions',
			marker=dict(size=marker_size,color='black', symbol='square')),
								row=2, col=1)





	fig.add_trace(go.Scatter(x=Disordered_Residues+1, y=Y_Protein,mode='markers',name='Protein Binding Regions',
			marker=dict(
					size=marker_size,
					color=Protein_DisordredScores,
					symbol='square',
					cmin=0.5,
					cmax=1.0,
					 #cmid=0.6,
					colorscale='Blues',
					showscale=False),
				meta=Protein_DisordredScores,
			 hovertemplate='Score: %{meta:.3f}'+'<br>Position: %{x}'),
								row=2, col=1)

	fig.add_trace(go.Scatter(x=Disordered_Residues+1, y=Y_DNA,mode='markers',name='DNA Binding Regions',
			marker=dict(
					size=marker_size,
					color=DNA_DisordredScores, #set color equal to a variable
					symbol='square',
					cmin=0.4,
					cmax=0.6,
					colorscale='Greens', # one of plotly colorscales
					showscale=False),
			 meta=DNA_DisordredScores,
			 hovertemplate='Score: %{meta:.3f}'+'<br>Position: %{x}'),
								row=2, col=1)

	fig.add_trace(go.Scatter(x=Disordered_Residues+1, y=Y_RNA,mode='markers',name='RNA Binding Regions',
			marker=dict(
					size=marker_size,
					color=RNA_DisordredScores, #set color equal to a variable
					symbol='square',
					cmin=0.5,
					cmax=0.6,
					colorscale='OrRd', # one of plotly colorscales
					showscale=False),
					meta=RNA_DisordredScores,
				 hovertemplate='Score: %{meta:.3f}'+'<br>Position: %{x}'),
								row=2, col=1)

	fig.add_trace(go.Scatter(x=Disordered_Residues+1, y=Y_Linker,mode='markers',name='Linker Regions',opacity=1,
			marker=dict(
					size=marker_size,
					color=Linker_DisordredScores, #set color equal to a variable
					symbol='square',
					cmin=1.1,
					cmax=1.5,
					colorscale='BuPu', # one of plotly colorscales
					showscale=False),
					meta=Linker_DisordredScores,
					hovertemplate='Score: %{meta:.3f}'+'<br>Position: %{x}'),
								row=2, col=1)

	tickvals = [2.5,3,3.5,4,5,6]
	ticktext = ['Linker Region', 'RNA Binding', 'DNA Binding', 'Protein Binding','Disorder','Sequence']

	fig.update_layout( {'height':600,
											'width':800,
											'title_text':"",
											# 'title_text':str(CAID_ID[j]),
											'xaxis':{'range': [0, len(Position)+1],'mirror':True},
											'xaxis2':{'range': [0, len(Position)+1],'mirror':True},
											'yaxis':{'range': [0, 1], 'dtick': 0.2,'mirror':True,'fixedrange':True,'title':'flDPnn diorder propensity'},
											'yaxis2':{'tickvals': tickvals, 'ticktext': ticktext,'mirror':True,'fixedrange':True},
											'paper_bgcolor':'rgba(255,255,255)',
											'plot_bgcolor':'rgba(0, 0, 0, 0)',
											'legend_title_text':'Click on legend to show/hide',
											'showlegend':True,
											 }
											)



	fig.update_xaxes(showline=True, linewidth=1, linecolor='black', gridcolor='black')
	fig.update_yaxes(showline=True, linewidth=1, linecolor='black', gridcolor='black')


	fig.show(config={'modeBarButtonsToRemove': ['select2d','lasso2d','zoomIn2d','zoomOut2d','hoverClosestCartesian','hoverCompareCartesian']})
	if warning_code in CAID_ID[j]:
		fname = str(CAID_ID[j]).replace(warning_code, "")
	else:
		fname = str(CAID_ID[j])
	fig.write_html(output_dir+"/"+fname+'.html',config={'modeBarButtonsToRemove': ['select2d','lasso2d','zoomIn2d','zoomOut2d','hoverClosestCartesian','hoverCompareCartesian']})
	


# In[39]:


def Xfunction(Input_List):
	b=0
	for b in range(0,len(CAID_ID),1):
		try:
			ordered_indexes= np.concatenate(np.argwhere(DisorderBinary_Prediction_byProteins[b]==0))
		except:
			ordered_indexes= list(np.argwhere(DisorderBinary_Prediction_byProteins[b]==0))
		target_list=np.asarray(Input_List[b]) # getting annotation list that we gonna modify
		target_list=target_list.astype(str)#conevrt that target list to an array of string
		insertion=np.full(len(ordered_indexes),'X') # creating an array of X to insert
		np.put(target_list, ordered_indexes,insertion) # insert X to un mapped ordered indexes in target list
		new_anno=target_list.tolist() # convert modified array back to a list
		Input_List[b]=new_anno # assign modified list back to the original list of lists


# In[40]:


Xfunction(ProteinPrediction_Score_byProteins)
Xfunction(DNAPrediction_Score_byProteins)
Xfunction(RNAPrediction_Score_byProteins)
Xfunction(LinkerPrediction_Score_byProteins)

Xfunction(ProteinBinary_Prediction_byProteins)
Xfunction(DNABinary_Prediction_byProteins)
Xfunction(RNABinary_Prediction_byProteins)
Xfunction(LinkerBinary_Prediction_byProteins)


# In[41]:
Lines=[]
Lines.append('*************************************************')
Lines.append('-------------File Format-------------------------')
Lines.append('line 1: >protein ID (when the pssm matrix cannot be calculated for an input sequence, '+warning_msg+' is added in front of the protein ID which means that the prediction is based on default/low quality PSSM which may lead to lower quality of the disorder predictions.)')
Lines.append('line 2: sequence')
Lines.append('line 3: binary disorder prediction (1 = disordered residue/amino acid; 0 ordered residue)')
Lines.append('line 4: disorder propensity (higher value denotes higher likelihood that a given residue is disordered)')
Lines.append('line 5: binary protein-binding prediction (1 = disordered protein-binding residue; 0 other disordered residue; X ordered residue)')
Lines.append('line 6: protein-binding propensity (higher value denotes higher likelihood that a given residue is disordered and binds proteins; X ordered residue)')
Lines.append('line 7: binary DNA-binding prediction (1 = disordered DNA-binding residue; 0 other disordered residue; X ordered residue)')
Lines.append('line 8: DNA-binding propensity (higher value denotes higher likelihood that a given residue is disordered and binds DNA; X ordered residue)')
Lines.append('line 9: binary RNA-binding prediction (1 = disordered RNA-binding residue; 0 other disordered residue; X ordered residue)')
Lines.append('line 10: RNA-binding propensity (higher value denotes higher likelihood that a given residue is disordered and binds RNA; X ordered residue)')
Lines.append('line 11: binary linker prediction (1 = disordered linker residue; 0 other disordered residue; X ordered residue)')
Lines.append('line 12: linker propensity (higher value denotes higher likelihood that a given residue is the disordered linker; X ordered residue)')
Lines.append('*************************************************')
b=0
for b in range(0,len(CAID_ID),1):
	raw_id = str(CAID_ID[b])
	if warning_code in raw_id:
		seqid = raw_id.replace(warning_code,warning_msg)
	else:
		seqid = raw_id
	Lines.append('>'+seqid)
	x = str(list(AA_Sequences_byProteins[b]))
	x = x[1:-1]
	# x = x.replace(",", "")
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)

	x = str(list(DisorderBinary_Prediction_byProteins[b]))
	x = x[1:-1]
	# x = x.replace(",", "")
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)
	x = str(list(DisorderPrediction_Score_byProteins[b]))
	x = x[1:-1]
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)

	x = str(list(ProteinBinary_Prediction_byProteins[b]))
	x = x[1:-1]
	# x = x.replace(",", "")
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)
	x = str(list(ProteinPrediction_Score_byProteins[b]))
	x = x[1:-1]
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)

	x = str(list(DNABinary_Prediction_byProteins[b]))
	x = x[1:-1]
	# x = x.replace(",", "")
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)
	x = str(list(DNAPrediction_Score_byProteins[b]))
	x = x[1:-1]
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)

	x = str(list(RNABinary_Prediction_byProteins[b]))
	x = x[1:-1]
	# x = x.replace(",", "")
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)
	x = str(list(RNAPrediction_Score_byProteins[b]))
	x = x[1:-1]
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)

	x = str(list(LinkerBinary_Prediction_byProteins[b]))
	x = x[1:-1]
	# x = x.replace(",", "")
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)
	x = str(list(LinkerPrediction_Score_byProteins[b]))
	x = x[1:-1]
	x = x.replace(" ", "")
	x = x.replace("'", "")
	Lines.append(x)


# In[42]:


with open(output_dir+'/function_results.txt', 'w') as f:
	for item in Lines:
		f.write("%s\n" % item)

with open(output_dir+'/results.csv', 'w') as f:
	for item in Lines:
		f.write("%s\n" % item)


# In[ ]:
