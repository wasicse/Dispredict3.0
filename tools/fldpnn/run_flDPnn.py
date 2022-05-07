#!/usr/bin/env python3


"""
Author: Sina Ghadermarzi
email: sina.ghadermarzi@gmail.com
date: 04/05/2020



a python wrapper for flDPnn tool: takes a single/multi sequence fasta file and saves csv result in the same folder
usage:
python run_flDPnn.py <input fasta path>

"""


import os
import sys
import os
from datetime import datetime
from random import choice
from string import digits
import pandas
import numpy
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
# import io

flDPnn_folder = this_script_directory = os.path.dirname(os.path.realpath(__file__))
warning_code = "[WARNING]"
warning_msg = "Warning: This prediction is based on default/low quality PSSM which may lead to lower quality of the disorder predictions."
nodisorder_msg = "Prediction for this sequence contains no disorder"

DISORDPBIND_PROT_THRESHOLD = 0.7120000000000001
FMORPHPRED_THRESHOLD = 0.3111
DISORDPBIND_DNA_THRESHOLD = 0.5329999999999999
DISORDPBIND_RNA_THRESHOLD = 0.214
DFLPRED_THRESHOLD = 0.22699999999999998

MIN_REGION_LENGTH = 4




def read_mainpred(res_file_address):
	pred_df = pandas.read_csv(res_file_address, sep = "\t", header=None, names=['Residue Number', 'Residue Type', "Predicted Score for Disorder", "Binary Prediction for Disorder"])
	return pred_df

def read_disordpbind(file_address):
	with open (file_address) as f:
		spl = f.read().split("\n")
		r_strs = spl[3].split(" ")
		if "" in r_strs:
			r_strs.remove("")
		rdp_r = [float(x) for x in r_strs]
		d_strs = spl[5].split(" ")
		if "" in d_strs:
			d_strs.remove("")
		rdp_d = [float(x) for x in d_strs]
		p_strs = spl[7].split(" ")
		if "" in p_strs:
			p_strs.remove("")
		rdp_p = [float(x) for x in p_strs]
	return rdp_r,rdp_d,rdp_p

def read_dflpred(file_address):
	with open (file_address) as f:
		spl = f.read().split("\n")
		dfl_strs = spl[2].split(" ")
		dfl = [float(x) for x in dfl_strs]
	return dfl

def read_fmorfpred(file_address):
	with open(file_address) as f:
		spl = f.read().split("\n")
		fmorf_strs = spl[2].split(" ")
		fmorf = [float(x) for x in fmorf_strs]
	return fmorf




def highlight(txt_str,col):
	out_str = r'<font color="'+col+'"><b>'+txt_str+'</b></font>'
	return out_str



def gen_taskid(prog_name):
	random_code= (''.join(choice(digits) for i in range(2)))
	time_string=  datetime.now().strftime("%Y%m%d%H%M%S")
	taskid = prog_name+"_"+time_string + random_code
	return taskid

def read_and_prepare_fasta(in_path,out_path):
	with open(in_path) as in_file:
		in_str = in_file.read().replace("\r","")
	spls = in_str.split(">")[1:]
	seqids = []
	seqs = []
	seq_idx = 0
	out_str = ""
	for spl in spls:
		lines = spl.split("\n")
		seqid  =lines[0]
		seq = "".join(lines[1:])
		pure_seq = ""
		valid_seq = (len(seq)<=5000)
		if valid_seq:
			for i in range(0, len(seq)):
				if seq[i] in "ABCDEFGHIKLMNPQRSTVWXYZ":
					pure_seq += seq[i]
				else:
					valid_seq = False
					print("Sequence "+seqid+" with an invalid character "+seq[i]+" at residue number "+str(i)+" was ignored!")
		else:
			print("Sequence "+seqid+" was ignored, because it is longer than 5000 residues!")
		if valid_seq:
			seqids.append(seqid)
			seqs.append(pure_seq)
			out_str+= ">"+str(seqid)+"\n"+pure_seq+"\n"
			seq_idx+=1
	with open(out_path,"w") as outf:
		outf.writelines(out_str)
	return seqids, seqs



def compose_header(seqid,regions):
	header_str =">"+seqid+"\n"
	i = 1
	for region in regions:
		header_str+= "Disordered region "+str(i)
		header_str+= " [from "+str(region["start"]+1)+" to "+str(region["end"])+"]"
		if (region["dna"] or region["rna"] or region["prot"] or region["linker"]):
			header_str +=" ("
			if (region["dna"]):
				header_str += "likely DNA-binding, "
			if (region["rna"]):
				header_str += "likely RNA-binding, "
			if (region["prot"]):
				header_str += "likely protein-binding, "
			if (region["linker"]):
				header_str += "likely linker, "
			header_str =header_str[:-2]+")"
		header_str += "\n"
		i+=1
	return header_str


def annotate_regions(result_df,rdp_r,rdp_d,rdp_p,dfl,fmorf):
	arr = result_df[result_df.columns[3]].astype(str).values
	# print(result_df["binary prediction of disorder"].str)
	arr = "".join(list(arr))
	on_region = False
	region_length = 0
	regions = []
	region_start = -1
	for i in range(len(arr)):
		if on_region:
			if arr[i] == "1":
				region_length += 1
			else:
				if region_length >= MIN_REGION_LENGTH:
					regions.append({"start": region_start, "end": region_start + region_length})
				region_start = -1
				region_length = 0
				on_region = False

		else:
			if arr[i] == "1":
				on_region = True
				region_start = i
				region_length = 1
	if on_region:
		if region_length >= MIN_REGION_LENGTH:
			regions.append({"start": region_start, "end": region_start + region_length})
	for reg in regions:
		s = reg["start"]
		e = reg["end"]
		if numpy.mean(dfl[s:e])>DFLPRED_THRESHOLD:
			reg["linker"] = True
			result_df.loc[s:e-1,"Linker"] = "L"
		else:
			reg["linker"] = False

		if numpy.mean(rdp_d[s:e]) > DISORDPBIND_DNA_THRESHOLD:
			reg["dna"] = True
			result_df.loc[s:e-1,"DNA-binding"] = "D"
		else:
			reg["dna"] = False

		if numpy.mean(rdp_r[s:e])>DISORDPBIND_RNA_THRESHOLD:
			reg["rna"] = True
			result_df.loc[s:e-1,"RNA-binding"] = "R"
		else:
			reg["rna"] = False

		if numpy.mean(rdp_p[s:e])>DISORDPBIND_PROT_THRESHOLD or numpy.mean(fmorf[s:e-1])>FMORPHPRED_THRESHOLD:
			reg["prot"] = True
			result_df.loc[s:e-1,"Protein-binding"] = "P"
		else:
			reg["prot"] = False
	return regions


# def table2htmltag(df,header):
# 	table_tag = r'''
# <table class="table table-condensed">
# <tr>
# <th colspan="8" class="seq" style="background-color:#FFA500">''' + header.replace("\n", "<br>") + r'''</th>
# </tr>
# <tr>
# <th>''' +	highlight(df.columns[0],"black")	+ r'''</th>
# <th>''' +	highlight(df.columns[1],"black")	+ r'''</th>
# <th>''' +	highlight(df.columns[2],"black")	+ r'''</th>
# <th>''' +	highlight(df.columns[3],"black")	+ r'''</th>
# <th>''' +	highlight(df.columns[4],"green")	+ r'''</th>
# </tr>'''
# 	for idx,row in df.iterrows():
# 		if row[df.columns[3]] == "1":
# 			table_tag += r'<tr>' + '\n'
# 			table_tag += r'<td>' +highlight(row[df.columns[0]],"black") + r'</td>' + '\n'
# 			table_tag += r'<td class="aa">' + highlight(row[df.columns[1]],"black") + r'</td>' + '\n'
# 			table_tag += r'<td>'+	highlight(row[df.columns[2]],"black")	+ r'</td>' + '\n'
# 			table_tag += r'<td>'+	highlight(row[df.columns[3]],"black")	+ r'</td>' + '\n'
# 			table_tag += r'<td>'+	highlight(row[df.columns[4]],"green")	+ r'</td>' + '\n'
# 			table_tag += r'</tr>' + '\n'
# 		else:
# 			table_tag += r'<tr>' + '\n'
# 			table_tag += r'<td>' + row[df.columns[0]] + r'</td>' + '\n'
# 			table_tag += r'<td class="aa">' + row[df.columns[1]] + r'</td>' + '\n'
# 			table_tag += r'<td>' + row[df.columns[2]] + r'</td>' + '\n'
# 			table_tag += r'<td>' + row[df.columns[3]] + r'</td>' + '\n'
# 			table_tag += r'<td>' + row[df.columns[4]] + r'</td>' + '\n'
# 			table_tag += r'</tr>' + '\n'
# 	table_tag += "</table>"
# 	return table_tag


def write_to_html(seqids,warnings,nodisorder_notice,table_df_list,out_path):
	html_begining = r'''
	<!DOCTYPE html>
	<html lang="en">
	<head>
	<meta charset="utf-8">
	<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
	<meta http-equiv="X-UA-Compatible" content="IE=Edge">
	<title>flDPnn Server - Results Page</title>
	</head>
	<body>
	<div class="container">
	<div class="row">
	<div class="col-lg-8 col-lg-offset-2">
	<h1>flDPnn results page</h1>
	</div>
	</div>
	<div class="row">
	<div class="col-lg-8 col-lg-offset-2">
	<div class="Predictions">
	<h2>Results</h2>
	<div class="table-responsive">
'''
	html_end  = r'''
	</div>
	</div>
	</div>
	<div class="row">
	<div class="col-lg-8 col-lg-offset-2">
	<footer>
	<h2>Visit biomine lab web page</h2>
	<a target="_blank" href="http://biomine.cs.vcu.edu">http://biomine.cs.vcu.edu</a>
	</footer>
	</div>
	</div>
	</div>
	</body>
	</html>
	'''
	html_middle =""
	for i in range(len(seqids)):
		html_middle += '<h2>>'+seqids[i]+"</h2>"
		if warnings[i]:
			html_middle += '<p style="color:red;">'+warning_msg+"</p>"
		if nodisorder_notice[i]:
			html_middle += '<p>'+nodisorder_msg+"</p>"
		html_middle += '''<object width="900px" height = "600px" data="'''+seqids[i]+".html"+r'''"></object>'''
	with open(out_path,"w") as outf:
		outf.writelines(html_begining)
		outf.writelines(html_middle)
		outf.writelines(html_end)















in_path = sys.argv[1]

taskid = gen_taskid("pyflDPnn_tmp")
# with open ("command_"+taskid+".txt","w") as commandfile:
# 	commandfile.writelines(in_path+"\n"+out_path)

in_path = os.path.abspath(in_path)
inout_dir = os.path.dirname(in_path)


cwd = os.getcwd()
os.chdir(flDPnn_folder)
os.mkdir(taskid)



# os.system("cp "+in_path+" " + flDPnn_folder + "/" + taskid + "/seqs.fasta")

seqids, seqs = read_and_prepare_fasta(in_path, flDPnn_folder + "/" + taskid + "/seqs.fasta")



os.system("./DisoComb.sh " + flDPnn_folder + "/" + taskid + "/seqs.fasta > /dev/null")

# os.system("cp "+flDPnn_folder+"/"+taskid+"/result.csv "+out_path)



# res_df_list = []
# header_list = []
# # txt_out_str = "Residue Number,Residue Type,Disorder Propensity, Binary Prediction of disorder"
# # for seqid in seqids:
# txt_out_str= ""

# warnings = [None]*len(seqids)
# nodis_flags = [None]*len(seqids)

# for i in range(len(seqids)):
# 	length = len(seqs[i])
# 	seqid = seqids[i]
# 	warnings[i] = os.path.exists(flDPnn_folder + "/" + taskid + "/pred_seqs/use_default_pssm_seq"+str(i))
# 	if warnings[i]:
# 		seqid+= warning_code

# 	result_df = read_mainpred(flDPnn_folder + "/" + taskid + "/pred_seqs/seq"+str(i)+".nn.pred")
# 	nodis_flags[i] = (result_df['Binary Prediction for Disorder'] == 0).all()
# 	rdp_r,rdp_d, rdp_p= read_disordpbind(flDPnn_folder + "/" + taskid + "/tmp_seqs/seq"+str(i)+".rdp")
# 	dfl = read_dflpred(flDPnn_folder + "/" + taskid + "/tmp_seqs/seq"+str(i)+".dfl")
# 	fmorf = read_fmorfpred(flDPnn_folder + "/" + taskid + "/tmp_seqs/seq"+str(i)+".fmorf")
# 	rdp_r = numpy.array(rdp_r)
# 	rdp_d = numpy.array(rdp_d)
# 	rdp_p = numpy.array(rdp_p)
# 	dfl = numpy.array(dfl)
# 	fmorf = numpy.array(fmorf)
# 	dna_col = numpy.array([""]*length)
# 	result_df["DNA-binding"] = dna_col
# 	rna_col = numpy.array([""] * length)
# 	result_df["RNA-binding"] = rna_col
# 	prot_col = numpy.array([""] * length)
# 	result_df["Protein-binding"] = prot_col
# 	linker_col = numpy.array([""] * length)
# 	result_df["Linker"] = linker_col
# 	# regions = annotate_regions(result_df,rdp_r,rdp_d,rdp_p,dfl,fmorf)
# 	# entry_header = compose_header(seqid,regions)
# 	result_df = result_df.applymap(str)
# 	# result_df["Predicted Function"] = result_df["DNA-binding"]+result_df["RNA-binding"]+result_df["Protein-binding"]+result_df["Linker"]
# 	result_df.drop(labels = ["DNA-binding","RNA-binding","Protein-binding","Linker"],axis = 1,inplace=True)
# 	result_df["rdp_r"] = rdp_r
# 	result_df["rdp_d"] = rdp_d
# 	result_df["rdp_p"] = rdp_p
# 	result_df["dfl"] = dfl
# 	result_df["fmorf"] = fmorf
# 	res_df_list.append(result_df)
# 	# header_list.append(entry_header)
# 	s = StringIO()
# 	result_df.to_csv(s,index= False)
# 	res_df_str = s.getvalue()
# 	txt_out_str += ">"+ seqid +"\n"+ res_df_str


# with open(inout_dir+"/results.csv","w") as txt_outf:
# 	txt_outf.writelines(txt_out_str)

# os.system("python3 "+flDPnn_folder+"/flDPnn_Function_workflow.py "+inout_dir+"/results.csv")
# write_to_html(seqids,warnings,nodis_flags, res_df_list,inout_dir+"/results.html")

os.system("rm -r "+flDPnn_folder+"/"+taskid)
#
os.chdir(cwd)
