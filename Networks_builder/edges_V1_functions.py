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
	print(source['pop_name'],target['pop_name'])
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


if __name__ == '__main__':
	dict_path="../Additional_data/dict_v1_nodes.json"
	net_layers,dataframes=add_nodes_V1_in_nrrd (dict_path,1) 
	df_connection_info="../Additional_data/l4_to_l4_connections_info.csv"
	dict_types=dict()
	dict_types["exc"]=["exc1","exc2","exc3","exc4","exc5","exc6","exc7"]
	dict_types["pvalb"]=["PV1","PV2"]
	dict_types["sst"]=["SST1","SST2","SST3"]
	dict_types["vip"]=["VIP1","VIP2","VIP3","VIP4"]
	dict_types["htr3a"]=["htr3a"]
	net=distance_edges_within(net_layers[2],df_connection_info,dict_types,n_synapses)
	#net.build() #all connections within l4

	






