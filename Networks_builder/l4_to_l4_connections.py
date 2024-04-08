import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime
from datetime import timedelta
from scipy.special import erf
from bmtk.builder import NetworkBuilder
from bmtk.builder.auxi.node_params import positions_columnar,positions_nrrd,xiter_random
import matplotlib.pyplot as plt
import h5py
import time
from bmtk.utils.create_environment import create_environment
from bmtk.analyzer.spike_trains import plot_raster
from bmtk.analyzer.compartment import plot_traces
import json
import os

def distance (x_pre,x_post,z_pre,z_post) :
	d=np.sqrt((x_pre-x_post)**2+(z_pre-z_post)**2)
	return(d)

def gaussian_function (d,amplitude,mu) :
	p=amplitude*np.exp(-d**2/(2*mu**2))
	return(p)

def exp_function (d,amplitude,mu):
	p=amplitude*np.exp(-d/mu)
	return(p)

def distance_proba (source,target,amplitude,mu,fun_type) :
	#print(source)
	print(source["pop_name"],target["pop_name"])
	d=distance(source['positions'][0],target['positions'][0],source['positions'][2],target['positions'][2])
	if fun_type=="gauss" :
		proba=gaussian_function(d,amplitude,mu)
	else :
		if fun_type=="exp" :
			proba=exp_function(d,amplitude,mu)
			print(proba)
			if proba>100 :
				proba=100
	p_connect=np.random.binomial(1,proba/100)
	if p_connect == 1:
		#print(source["node_id"],target["node_id"],proba)
		return(1)
	else :
		return(0)

def billeh_gaussian_function (d,amplitude,mu) : 

	p=amplitude*np.exp(-d**2/(2*mu**2))

	return(p)


def distance_proba_billeh(source,target,amplitude,mu) : 
	print(source["pop_name"],target["pop_name"])
	d=distance(source['positions'][0],target['positions'][0],source['positions'][2],target['positions'][2])

	proba=billeh_gaussian_function(d,amplitude,mu)

	p_connect=np.random.binomial(1,proba/100)

	if p_connect == 1: 

		#print(source["node_id"],target["node_id"],proba)

		return(1)

	else : 

		return(0)

startt=time.time()
net=NetworkBuilder('cortex_l4_factor_0.1')

net.import_nodes(nodes_file_name='../Networks/nodes/cortex_l4_0.1_nodes.h5',node_types_file_name='../Networks/nodes/cortex_l4_0.1_node_types.csv')

l4_info_connections=pd.read_csv("../Additional_data/df_l4.csv")
df_htr3a=pd.read_csv("../Additional_data/htr3a_l4.csv",sep=";")

types=["exc","pvalb","sst","vip","htr3a"]
subtypes=[["exc1","exc2","exc3","exc4","exc5","exc6","exc7"],["PV1","PV2"],["SST1","SST2","SST3"],["VIP1","VIP2","VIP3","VIP4"],["htr3a"]]

for i in np.arange(l4_info_connections.shape[0]) : 
	pre_type=l4_info_connections.loc[i]["pre"]
	post_type=l4_info_connections.loc[i]["post"] 
	index_pre=types.index(pre_type)
	index_post=types.index(post_type)
	if l4_info_connections.loc[i]["gaussian_fun"]<l4_info_connections.loc[i]["exponential_fun"]:
		pmax=round(l4_info_connections.loc[i]["gaussian_pmax"]*100,2)
		sigma=round(l4_info_connections.loc[i]["gaussian_mu"]*1e6,2)
		function="gauss"
	else:  
		pmax=round(l4_info_connections.loc[i]["exponential_pmax"]*100,2)
		sigma=round(l4_info_connections.loc[i]["exponential_mu"]*1e6,2)
		function="exp"
	if pre_type=="exc" :
		if post_type=="exc":
			synaptic_type='static_ExchToExc.json'
			weight=10
		else : 
			synaptic_type='static_ExcToInh.json'
			weight=10
	else : 
		if post_type=="exc":
			synaptic_type='static_InhToExc.json'
			weight=-10
		else : 
			synaptic_type='static_InhToInh.json'
			weight=-10
	for pre in subtypes[index_pre] : 
		for post in subtypes[index_post] : 
			print(pre,post,pmax,sigma,function,synaptic_type,weight)
			start=time.time()
			net.add_edges(
				source={'pop_name':pre},
				target={'pop_name':post},
				connection_rule=distance_proba,
				connection_params={'amplitude':pmax,'mu':sigma,'fun_type':function},
        			syn_weight=weight,
				delay= np.random.uniform(1,3),
				dynamics_params=synaptic_type,
				model_template='static_synapse'
			)
			end=time.time()
			print(pre,post,end-start)

for i in np.arange(df_htr3a.shape[0]) : 
	pre_type=df_htr3a.loc[i]["pre"]
	post_type=df_htr3a.loc[i]["post"]
	pmax=round(df_htr3a.loc[i]["pmax"]*100,2)
	sigma=round(df_htr3a.loc[i]["mu"],2)
	index_pre=types.index(pre_type)
	index_post=types.index(post_type)
	if pre_type =="exc" : 
		synaptic_type='static_ExcToInh.json'
		weight=10
	else : 
		if post_type=="exc" : 
			synaptic_type='static_InhToExc.json'
			weight=-10
		else : 
			synaptic_type='static_InhToInh.json'
			weight=-10
	for pre in subtypes[index_pre] : 
		for post in subtypes[index_post] : 
			#print(pre,post,pmax,sigma,synaptic_type,weight)
			start=time.time()
			net.add_edges(
				source={'pop_name':pre},
				target={'pop_name':post},
				connection_rule=distance_proba_billeh,
				connection_params={'amplitude':pmax,'mu':sigma},
				syn_weight=weight,
				delay= np.random.uniform(1,3),
				dynamics_params=synaptic_type,
				model_template='static_synapse'
			)
			end=time.time()
			print(pre,post,end-start)
			

net.build()
net.save_edges(output_dir='../Networks/l4_to_l4')

endd=time.time()
print(endd-startt)

