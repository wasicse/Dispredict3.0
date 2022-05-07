import numpy as np
import os
import sys
import tensorflow as tf
import keras
from keras import backend as K
from keras.models import load_model
thres = 0.505
modelname = "layer64_8.2.h5"
def readfea(filename):
	data=np.loadtxt(filename)
	plen=data.shape[0]
	features=np.zeros([data.shape[0],data.shape[1]+3])
	for i in range(plen):
		tt = data[i].tolist()
		tt.append(plen)
		tt.append(i)
		tt.append(min(i,plen-i-1))
		features[i] = tt
	return features

feas=readfea(sys.argv[1])
model=load_model(modelname)
preds = model.predict(feas,verbose=0)
bins = [item>=thres for item in preds]
for p,b in zip(preds,bins):
	print("%1.3f\t%d" % (p,b))
