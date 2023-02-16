# Author : Md Wasi Ul Kabir  

import torch
import esm
import pathlib 
import joblib
import numpy as np
import warnings
from Bio import SeqIO
import gdown
import subprocess
from optparse import OptionParser
import os
import random
from pathlib import Path

warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  
np.set_printoptions(precision=3)

# Set a seed value: 
seed_value= 2515  
os.environ['PYTHONHASHSEED']=str(seed_value) 
random.seed(seed_value) 
np.random.seed(seed_value) 

def loadModels():
    print("Loading models...")
    output = parent_path+"/models/model.pkl"
    path = pathlib.Path(output)
    if not path.is_file():
        
        url = "https://www.cs.uno.edu/~mkabir3/Dispredict3.0/models/model.pkl"
        gdown.download(url=url, output=output, quiet=False, fuzzy=True)

    output = parent_path+"/models/pca.pkl"
    path = pathlib.Path(output)
    if not path.is_file():
        url = "https://www.cs.uno.edu/~mkabir3/Dispredict3.0/models/pca.pkl"
        gdown.download(url=url, output=output, quiet=False, fuzzy=True)

    output = parent_path+"/models/scaler.pkl"
    path = pathlib.Path(output)
    if not path.is_file():
        url = "https://www.cs.uno.edu/~mkabir3/Dispredict3.0/models/scaler.pkl"
        gdown.download(url=url, output=output, quiet=False, fuzzy=True)



    output = parent_path+"/tools/fldpnn/programs/blast-2.2.24/db/swissprot.psq"
    path = pathlib.Path(output)
    if not path.is_file():
        url = "https://www.cs.uno.edu/~mkabir3/Dispredict3.0/db/swissprot.psq"
        gdown.download(url=url, output=output, quiet=False, fuzzy=True)

    output = parent_path+"/tools/fldpnn/programs/blast-2.2.24/db/swissprot.phr"
    path = pathlib.Path(output)
    if not path.is_file():
        url = "https://www.cs.uno.edu/~mkabir3/Dispredict3.0/db/swissprot.phr"
        gdown.download(url=url, output=output, quiet=False, fuzzy=True)


def dispredict(fasta_filepath,output_path):
    
    path = pathlib.Path(output_path+fasta_filepath.split("/")[-1].split(".")[0]+"_disPred.txt")
    if path.is_file():
            print("Removing old output files...") 
            bashCommand="rm -rf "+output_path+fasta_filepath.split("/")[-1].split(".")[0]+"_disPred.txt"
            output = subprocess.check_output(["bash","-c", bashCommand])
            print(output.decode('utf-8')) 

            bashCommand="rm -rf "+parent_path+"/tools/fldpnn/output/*"
            output = subprocess.check_output(["bash","-c", bashCommand])
            print(output.decode('utf-8')) 

    print("Extracting features from fldpnn...") 
    bashCommand="python "+parent_path+"/tools/fldpnn/run_flDPnn.py "+fasta_filepath
    output = subprocess.check_output(["bash","-c", bashCommand])
    print(output.decode('utf-8')) 

    print("Loading ESM-1b model...")
    # model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")

    # Load ESM-1b model
    
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    batch_converter = alphabet.get_batch_converter()

    # %%
    print("Dispredict3.0 prediction started...")
    threshold=0.382
    # flagpid=0
    for record in SeqIO.parse(fasta_filepath, "fasta"):
        print(record.id)
        # pid=record.id.split("|")[1].()
        pid=record.id.strip()
        fasta=record.seq
        # print(pid)
        # print(fasta)

        print("Sequence length: ", len(fasta))

        if(len(fasta)>1022):
                n = 1022
                chunks = [fasta[i:i+n] for i in range(0, len(fasta), n)]
                # print("No. of Chunks: ",len(chunks))

                flagc=0
                for chunk in chunks:
                    # print("chunk length: ",len(chunk))
                    # print(chunk)

                    data = [ (pid, chunk)]
                    batch_labels, batch_strs, batch_tokens = batch_converter(data)
                    flagl=0
                    # print("Layer", layern)
                    layern=list(range(34)) 
                    with torch.no_grad():
                        results = model(batch_tokens, repr_layers=layern, return_contacts=True)
                    # Extract per-residue representations (on CPU)
                    for layern in range(34): #
                                                                        
                        token_representations = results["representations"][layern]
                        toekn=token_representations[0, 1 : len(chunk) + 1].numpy()
                        toekn=np.hstack((toekn,np.mean(toekn, axis=1).reshape(-1,1)))
                        if(flagl==0):
                            np_featurel=toekn
                            flagl=1
                        else:
                            np_featurel=np.hstack((np_featurel,toekn))                       
        

                    if(flagc==0):
                        np_featurell=np_featurel
                        flagc=1
                    else:
                        np_featurell=np.vstack((np_featurell,np_featurel))
                    # print(np_featurell.shape)

        else:

            data = [ (pid, fasta)]
            batch_labels, batch_strs, batch_tokens = batch_converter(data)
            flagl=0
            layern=list(range(34))    
            # print("Layer", layern)
            # Extract per-residue representations (on CPU)
            with torch.no_grad():
                results = model(batch_tokens, repr_layers=layern, return_contacts=True)
                        
            for layern in range(34): #               
                
                token_representations = results["representations"][layern]
                toekn=token_representations[0, 1 : len(fasta) + 1].numpy()
                toekn=np.hstack((toekn,np.mean(toekn, axis=1).reshape(-1,1)))

                if(flagl==0):
                    np_featurell=toekn
                    flagl=1
                else:
                    np_featurell=np.hstack((np_featurell,toekn))                  
        
        np_fldscore=np.loadtxt(parent_path+"/tools/fldpnn/output/"+pid+".ttscore")
        np_proba_pred=np.loadtxt(parent_path+"/tools/fldpnn/output/"+pid+".ttpreds")
        np_index=np.loadtxt(parent_path+"/tools/fldpnn/output/"+pid+".ttindex",dtype='object')


        np_fld=np.hstack((np_fldscore,np_proba_pred))       

        scaler= joblib.load(parent_path+"/models/scaler.pkl")
        np_featurell = scaler.transform(np_featurell)

        ipca= joblib.load(parent_path+"/models/pca.pkl")

        np_featurell = ipca.transform(np_featurell)

        np_allfeat=np.hstack((np_fld, np_featurell))   
       

        saved_model = joblib.load(parent_path+"/models/model.pkl")
        # print(np_allfeat.shape)
        proba = saved_model.predict_proba(np_allfeat)

        pred = (proba[:,1] >= threshold).astype(np.int)

        result=np.hstack((np_index,np.round(proba[:,1], 3).reshape(-1,1) ,pred.reshape(-1,1))) 


        

        num_ones = (pred == 1).sum()
        fd_proba=num_ones/pred.shape[0]
        if(fd_proba > 0.95):
                # print("Fully disorder proteins") 
                fd_label=1               
        else:
                # print("Not Fully disorder proteins") 
                fd_label=0 

        with open(output_path+fasta_filepath.split("/")[-1].split(".")[0]+"_disPred.txt", "ab") as f:
            f.write((">"+pid+"\n").encode())
            fmt = '%s', '%s', '%1.3f', '%s'
            np.savetxt(f, result, delimiter='\t',fmt=fmt) 
            
        with open(output_path+fasta_filepath.split("/")[-1].split(".")[0]+"_fullydisPred.txt", "ab") as f:
            f.write((">"+pid+"\n").encode())
            fmt = '%1.3f', '%s'
            np.savetxt(f, np.array([fd_proba, fd_label ]).reshape(1,2), delimiter='\t',fmt=fmt) 


    bashCommand="rm -rf "+parent_path+"tools/fldpnn/output/*"
    output = subprocess.check_output(["bash","-c", bashCommand])
    print(output.decode('utf-8')) 
    print("Dispredict3.0 prediction end...")
if __name__ == '__main__':
    
    parent_path = str(Path(__file__).resolve().parents[1])
    # print("Parent Dir",parent_path)
    
    parser = OptionParser()
    parser.add_option("-f", "--fasta_filepath", dest="fasta_filepath", help="Path to input fasta.", default=parent_path+'/example/sample.fasta')
    parser.add_option("-o", "--output_path", dest="output_path", help="Path to output.", default=parent_path+'/output/')

    (options, args) = parser.parse_args()

    print("Dataset Path:",options.fasta_filepath)
    print("Output Path:",options.output_path)

    workspace=options.output_path
    pathlib.Path(workspace).mkdir(parents=True, exist_ok=True) 

    workspace=parent_path+"/models"
    pathlib.Path(workspace).mkdir(parents=True, exist_ok=True)

    loadModels()
    
    dispredict(options.fasta_filepath,options.output_path)
    

