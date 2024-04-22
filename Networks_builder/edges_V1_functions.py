#!/usr/bin/env python3
from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder
import plotly.express as px
import math
import numpy as np
import json
from position_building_V1 import get_population_loc_df,create_layer_population_loc,add_nodes_V1_in_nrrd

def distance (x_pre,x_post,z_pre,z_post) :
	d=np.sqrt((x_pre-x_post)**2+(z_pre-z_post)**2)
	return(d)


def distance_connection(source,target,amplitude,mu,function_type,n_synapses) : #only works for now for nodes h5 with positions (instead of x,y,z) 
	
	#print(source['location'],source['pop_name'],source['node_id'],target["location"],target['pop_name'],target['node_id'])
	d=distance(source['positions'][0],target['positions'][0],source['positions'][2],target['positions'][2])
	if function_type == "gaussian" : 
		proba=amplitude*np.exp(-d**2/(2*mu**2))
	else : 
		if function_type=="exponential" : 
			proba=amplitude*np.exp(-d/mu)
		else : 
			proba=amplitude*np.exp(-d**2/(2*mu**2)) #billeh function
	p_connect=np.random.binomial(1,proba/100)
	if p_connect == 1:
		return(n_synapses)
	else : 
		return(0)

def distance_edges_within (net,df_connection_info,dict_types,n_synapses) : #faire une fonction qui choisit le nomber de n ou le mettre dans le dict
	for i in np.arange(df_connection_info.shape[0]) : 
		pre_type=df_connection_info.loc[i]["pre"]
		post_type=df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[pre_type] : 
			for post_subtype in dict_types[post_type] : 
				net.add_edges(
					source={'pop_name':pre_subtype},
					target={'pop_name':post_subtype},
					connection_rule=distance_connection,
					connection_params={'amplitude':df_connection_info.loc[i]["pmax"],'mu':df_connection_info.loc[i]["sigma"],'function_type':df_connection_info.loc[i]["rule"],'n_synapses':n_synapses},
        				syn_weight=df_connection_info.loc[i]["weight"],
					delay= np.random.uniform(1,3),
					dynamics_params=df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
	return(net)
	
	
def distance_edges_between (net_pre,net_post,df_connection_info_path,dict_types,n_synapses) : #faire une fonction qui choisit le nomber de n ou le mettre dans le dict (arugment devrait être un interval)
	#trouver une façon de prendre seulement les lignes du df_connections selon les layer de net_pre et net_post : avec getattr(net_pre._name) et split pour obtenir l+nber layer
	df=pd.read_csv(df_connection_info_path)
	name_pre=net_pre.name
	split_pre=name_pre.split("_")
	layer_pre="l"+str(split_pre[1])
	name_post=net_post.name
	split_post=name_post.split("_")
	layer_post="l"+str(split_post[1])
	sub_df_connection_info=df.loc[(df["pre_layer"]==layer_pre)&(df["post_layer"]==layer_post)].dropna()
	sub_df_connection_info=sub_df_connection_info.drop(["Unnamed: 0"],axis=1)
	for i in list(sub_df_connection_info.index) : 
		pre_type=sub_df_connection_info.loc[i]["pre"]
		post_type=sub_df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[layer_pre][pre_type] : 
			for post_subtype in dict_types[layer_post][post_type] : 
				net_pre.add_edges(
					source=net_pre.nodes(pop_name=pre_subtype),
					target=net_post.nodes(pop_name=post_subtype),
					connection_rule=distance_connection,
					connection_params={'amplitude':sub_df_connection_info.loc[i]["pmax"],'mu':sub_df_connection_info.loc[i]["sigma"],'function_type':sub_df_connection_info.loc[i]["rule"],'n_synapses':n_synapses},
        				syn_weight=sub_df_connection_info.loc[i]["weight"],
					delay= np.random.uniform(1,3),
					dynamics_params=sub_df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
	return(net_pre)
	
	

"""
if __name__ == '__main__':
	dict_path="../Additional_data/dict_v1_nodes.json"
	net_layers,dataframes=add_nodes_V1_in_nrrd (dict_path,1) 
	df_connection_info_path="../Additional_data/connection_infos/from_l1_connections_info.csv"
	dict_types=dict()
	#
	dict_types_pre["exc"]=["exc1","exc2","exc3","exc4","exc5","exc6","exc7"]
	dict_types_pre["pvalb"]=["PV1","PV2"]
	dict_types_pre["sst"]=["SST1","SST2","SST3"]
	dict_types_pre["vip"]=["VIP1","VIP2","VIP3","VIP4"]
	dict_types_pre["htr3a"]=["htr3a"]
	#
	dict_l1=dict()
	dict_l1["htr3a"]=["htr3a1","htr3a2"]
	dict_l23=dict()
	dict_l23["exc"]=["exc"]
	dict_l23["pv"]=["PV1","PV2","PV3","PV4"]
	dict_l23["sst"]=["SST1","SST2","SST3"]
	dict_l23["vip"]=["VIP1","VIP2","VIP3","VIP4","VIP5"]
	dict_l23["htr3a"]=["htr3a1","htr3a2"]
	dict_l4=dict()
	dict_l4["exc"]=["exc1","exc2","exc3","exc4","exc5","exc6","exc7"]
	dict_l4["pv"]=["PV1","PV2"]
	dict_l4["sst"]=["SST1","SST2","SST3"]
	dict_l4["vip"]=["VIP1","VIP2","VIP3","VIP4"]
	dict_l4["htr3a"]=["htr3a"]
	dict_l5=dict()
	dict_l5["cf"]=["CF"]
	dict_l5["it"]=["IT1","IT2"]
	dict_l5["pv"]=["PV"]
	dict_l5["sst"]=["SST1","SST2","SST3","SST4"]
	dict_l5["vip"]=["VIP"]
	dict_l5["htr3a"]=["htr3a"]
	dict_l6=dict()
	dict_l6["exc"]=["exc1","exc2","exc3"]
	dict_l6["pv"]=["PV1","PV2","PV3","PV4"]
	dict_l6["sst"]=["SST1","SST2","SST3"]
	dict_l6["vip"]=["VIP"]
	dict_l6["htr3a"]=["htr3a"]
	dict_types={"l1":dict_l1,"l23":dict_l23,"l4":dict_l4,"l5":dict_l5,"l6":dict_l6}
	
	n_synapses=10
	#net=distance_edges_within(net_layers[2],df_connection_info,dict_types,n_synapses)
	net=distance_edges_between(net_layers[0],net_layers[1],df_connection_info_path,dict_types,n_synapses) #changer argument
	#net.build() #all connections within l4
"""
	

	






