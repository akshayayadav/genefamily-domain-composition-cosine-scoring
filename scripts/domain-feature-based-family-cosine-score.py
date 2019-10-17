#!/usr/bin/python3

#Author: Akshay Yadav
#Version: 1.0.0

import re
import os
import operator
import sys
import argparse
import subprocess
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np

parser = argparse.ArgumentParser(description="Tool for calculating domain feature based cosine score for given set of gene families")
parser._optionals.title="Arguments"

parser.add_argument('--pfamscanout_dir', help="Location of the pfamscan output files", required=True, dest="pfamscan_directory")
parser.add_argument('--output', help="Location of the output file", required=True, dest="output_file_name")
args = parser.parse_args()

###############################################################################################################################################
def read_pfamscanout_file(pfamscanout_fileName, pfamscanout_dirName):
	pfamscanout_file=open(pfamscanout_dirName+"/"+pfamscanout_fileName,"r")
	seqid_startcoord_domain_score_dict={}
	for line in pfamscanout_file:
		line=line.strip()
		if(re.match(r'^\#',line) or re.match(r'^$',line)):
			continue;
		linearr=re.split(r'\s+',line)
		seqid=linearr[0]
		start_coord=int(linearr[3])
		bit_score=linearr[11]
		domain_name=linearr[6]
		clan_id=linearr[14]
		#if not (re.match(r'^No\_clan$',clan_id)):
		#	domain_name=clan_id

		if(seqid in seqid_startcoord_domain_score_dict):
			seqid_startcoord_domain_score_dict[seqid][start_coord]={}
			seqid_startcoord_domain_score_dict[seqid][start_coord][domain_name]=bit_score
		else:
			seqid_startcoord_domain_score_dict[seqid]={}
			seqid_startcoord_domain_score_dict[seqid][start_coord]={}
			seqid_startcoord_domain_score_dict[seqid][start_coord][domain_name]=bit_score
	
	pfamscanout_file.close()
	
	return(seqid_startcoord_domain_score_dict)

def process_sequence_pairs(seqid_startcoord_domain_score_dict, pfamscanout_fileName, output_file):
	seqid_arr = list(seqid_startcoord_domain_score_dict.keys())
	cosine_score_arr=list()
	famsize=len(seqid_arr)
	if(famsize<2):
		return(0)
	for i in range(0,famsize):
		for j in range(i+1,famsize):
			startcoord_domain_score_dict_seq1=seqid_startcoord_domain_score_dict[seqid_arr[i]]
			startcoord_domain_score_dict_seq2=seqid_startcoord_domain_score_dict[seqid_arr[j]]
			cosine_score=calculate_domain_feature_cosine_score(startcoord_domain_score_dict_seq1, startcoord_domain_score_dict_seq2)
			cosine_score_arr.append(cosine_score)
			#print '{0}\t{1}\t{2}'.format (seqid_arr[i], seqid_arr[j],cosine_score)
	mean_cosine_score=calculate_mean_cosine_score(cosine_score_arr)
	output_file.write(os.path.splitext(os.path.basename(pfamscanout_fileName))[0]+" "+str(famsize)+" "+ str(mean_cosine_score[0])+"\n")

def calculate_domain_feature_cosine_score(startcoord_domain_score_dict_seq1, startcoord_domain_score_dict_seq2):
	startcoord_domain_score_dict_seq1_sorted=sorted(startcoord_domain_score_dict_seq1.items(), key=operator.itemgetter(0))
	startcoord_domain_score_dict_seq2_sorted=sorted(startcoord_domain_score_dict_seq2.items(), key=operator.itemgetter(0))

	
	seq1_domain_id_domain_score_dict=get_domain_features(startcoord_domain_score_dict_seq1_sorted)
	seq2_domain_id_domain_score_dict=get_domain_features(startcoord_domain_score_dict_seq2_sorted)

	feature_id_arr= get_feature_ids(seq1_domain_id_domain_score_dict, seq2_domain_id_domain_score_dict)

	cosine_score=get_cosine_score_from_feature_vector(seq1_domain_id_domain_score_dict, seq2_domain_id_domain_score_dict, feature_id_arr)

	return(cosine_score)
	

def get_domain_features(startcoord_domain_score_dict_seq_sorted):
	duplicate_counter=2
	domain_id_domain_score_dict={}
	feature_dict={}
	for entry in startcoord_domain_score_dict_seq_sorted:
		domain_id=list(entry[1].keys())[0]
		domain_score=list(entry[1].values())[0]
		if(domain_id in feature_dict):
			domain_id_domain_score_dict[domain_id+"-"+str(duplicate_counter)]=domain_score
			duplicate_counter+=1
		else:
			feature_dict[domain_id]=1
			domain_id_domain_score_dict[domain_id]=domain_score
	
	return(domain_id_domain_score_dict)

def get_feature_ids(seq1_domain_id_domain_score_dict, seq2_domain_id_domain_score_dict):
	feature_id_arr=list(seq1_domain_id_domain_score_dict.keys()) + list(seq2_domain_id_domain_score_dict.keys())
	feature_id_arr=list(set(feature_id_arr))
	return(feature_id_arr)

def get_cosine_score_from_feature_vector(seq1_domain_id_domain_score_dict, seq2_domain_id_domain_score_dict, feature_id_arr):
	
	seq1_feature_arr=get_feature_vector(seq1_domain_id_domain_score_dict, feature_id_arr)
	seq2_feature_arr=get_feature_vector(seq2_domain_id_domain_score_dict, feature_id_arr)

	seq1_feature_arr = np.array(seq1_feature_arr).reshape(1, -1)
	seq2_feature_arr = np.array(seq2_feature_arr).reshape(1, -1)
	cosine_score=cosine_similarity(seq1_feature_arr, seq2_feature_arr)
	cosine_score = cosine_score[0][0]
	return(cosine_score)
	
def get_feature_vector(seq_domain_id_domain_score_dict, feature_id_arr):
	seq_feature_arr=list()
	
	for feature_id in feature_id_arr:
		if(feature_id in seq_domain_id_domain_score_dict):
			seq_feature_arr.append(seq_domain_id_domain_score_dict[feature_id])
		else:
			seq_feature_arr.append(0)
	return(seq_feature_arr)

def calculate_mean_cosine_score(cosine_score_arr):
	cosine_score_arr=np.array(cosine_score_arr).reshape(1, -1)
	mean_cosine_score=np.mean(cosine_score_arr, axis=1)
	return(mean_cosine_score)

def run_workflow(output_fileName, pfamscanout_dirName):

	if not os.path.exists(os.path.dirname(output_fileName)):
    		os.makedirs(os.path.dirname(output_fileName))

	output_file=open(output_fileName, "w")
	for pfamscanout_fileName in os.listdir(pfamscanout_dirName):
		seqid_startcoord_domain_score_dict = read_pfamscanout_file(pfamscanout_fileName, pfamscanout_dirName)
		process_sequence_pairs(seqid_startcoord_domain_score_dict, pfamscanout_fileName, output_file)
	output_file.close()
		
#################################################################################################################################################
pfamscanout_dirName = args.pfamscan_directory 
output_fileName=args.output_file_name
run_workflow(output_fileName, pfamscanout_dirName)
