#!/usr/bin/env python3
from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder
import math
import numpy as np
import json
from position_building_V1 import get_population_loc_df,create_layer_population_loc,add_nodes_V1_in_nrrd
from edges_V1_functions import distance_edges_between,lgn_to_v1_layer,distance_edges_between_with_custom_delay,distance_edges_within,distance_edges_within_with_custom_delay,distance_edges_within_orientation,distance_edges_between_orientation,distance_orientation_custom_delay_between,distance_orientation_custom_delay_within 

def velocity_extraction(velocity_path,layer_pre,layer_post) :
	df=pd.read_csv(velocity_path)
	velocity=list(df.loc[df["Unnamed: 0"]==layer_pre][layer_post])[0]
	return(velocity)

def build_network_within (path_h5,path_csv,df_connection_info_path,dict_types,n_synapses,output_dir):
	path_split=path_h5.split("/")
	path_split_bis=path_split[-1].split("_nodes")
	net=NetworkBuilder(path_split_bis[0])
	net.import_nodes(nodes_file_name=path_h5,node_types_file_name=path_csv)
	network=distance_edges_within(net,df_connection_info_path,dict_types,n_synapses)
	network.build()
	network.save_edges(output_dir=output_dir)
    
def build_network_between (dict_node_path,factor,df_connection_info_path,dict_types,n_synapses,id_pre_layer,id_post_layer,output_dir) : 
	net_layers,dataframes= add_nodes_V1_in_nrrd(dict_node_path,factor)
	network=distance_edges_between(net_layers[id_pre_layer],net_layers[id_post_layer],df_connection_info_path,dict_types,n_synapses)
	network.build()
	network.save(output_dir) 
	
	
	
def build_network_between_with_existing_nodes (path_pre_h5,path_pre_csv,path_post_h5,path_post_csv,df_connection_info_path,dict_types,n_synapses,output_dir,custom_delay,path_velocity,orientation,slope) :
	path_pre_split=path_pre_h5.split("/")
	path_pre_split_bis=path_pre_split[-1].split("_nodes")
	net_pre=NetworkBuilder(path_pre_split_bis[0])
	path_post_split=path_post_h5.split("/")
	path_post_split_bis=path_post_split[-1].split("_nodes")
	net_post=NetworkBuilder(path_post_split_bis[0])
	net_pre.import_nodes(nodes_file_name=path_pre_h5,node_types_file_name=path_pre_csv)
	net_post.import_nodes(nodes_file_name=path_post_h5,node_types_file_name=path_post_csv)
	if custom_delay ==True :
		if path_pre_h5==path_post_h5 : # Inside a same population
			split_name=path_post_split_bis[0].split("_")
			if orientation == True : 
				network=distance_orientation_custom_delay_within(net_pre,df_connection_info_path,dict_types,n_synapses,slope,path_velocity,split_name[1])
			else :
				network=distance_edges_within_with_custom_delay(net_pre,df_connection_info_path,dict_types,n_synapses,path_velocity,split_name[1])
			
		else :
			if orientation == True : 
				network=distance_orientation_custom_delay_between(net_pre,net_post,df_connection_info_path,dict_types,n_synapses,slope,path_velocity)
			else : 
				network=distance_edges_between_with_custom_delay (net_pre,net_post,df_connection_info_path,dict_types,n_synapses,path_velocity)
			
	else: 
		if path_pre_h5==path_post_h5 : 
			if orientation == True : 
				network=distance_edges_within_orientation (net_pre,df_connection_info_path,dict_types,n_synapses,slope)
			else : 
				network=distance_edges_within(net_pre,df_connection_info_path,dict_types,n_synapses)
		else:
			if orientation == True : 
				network=distance_edges_between_orientation(net_pre,net_post,df_connection_info_path,dict_types,n_synapses,slope)
			else : 
				network=distance_edges_between(net_pre,net_post,df_connection_info_path,dict_types,n_synapses)
	network.build()
	network.save_edges(output_dir=output_dir)  

def build_network_between_lgn_v1 (path_pre_h5,path_pre_csv,path_post_h5,path_post_csv,lgn_to_l4_dict,layer_types,field_size,output_dir) :
	path_pre_split=path_pre_h5.split("/")
	path_pre_split_bis=path_pre_split[-1].split("_nodes")
	net_pre=NetworkBuilder(path_pre_split_bis[0])
	path_post_split=path_post_h5.split("/")
	path_post_split_bis=path_post_split[-1].split("_nodes")
	net_post=NetworkBuilder(path_post_split_bis[0])
	net_pre.import_nodes(nodes_file_name=path_pre_h5,node_types_file_name=path_pre_csv)
	net_post.import_nodes(nodes_file_name=path_post_h5,node_types_file_name=path_post_csv)
	network=lgn_to_v1_layer(net_pre,net_post,lgn_to_l4_dict,layer_types,field_size)
	network.build()
	network.save_edges(output_dir=output_dir) 	


"""
if __name__ == '__main__':
	lgn_to_l4_path="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/Additional_data/lgn_to_l4_dict.json"
	lgn_to_l4_dict=json.load(open(lgn_to_l4_path,"r"))
	layer_types=["exc1","exc2","exc3","exc4","exc5","exc6","exc7","PV1","PV2","SST1","SST2","SST3","VIP1","VIP2","VIP3","VIP4","htr3a"]
	field_size=(240.0, 120.0)
	output_dir="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/stockage_linux/lgn/100p/connections"
	pre_path_h5="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/stockage_linux/lgn/100p/lgn_nodes.h5"
	pre_path_csv="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/stockage_linux/lgn/100p/lgn_node_types.csv"
	post_path_h5="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/stockage_linux/network_v1_01/nodes/layer_4_factor_1_nodes.h5"
	post_path_csv="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/stockage_linux/network_v1_01/nodes/layer_4_factor_1_node_types.csv"
	build_network_between_lgn_v1(pre_path_h5,pre_path_csv,post_path_h5,post_path_csv,lgn_to_l4_dict,layer_types,field_size,output_dir)
"""		

"""
if __name__ == '__main__':
	dict_node_path="../Additional_data/dict_v1_nodes_bis.json"
	factor=0.1
	df_connection_info_path="../Additional_data/connection_infos/from_l4_connections_info.csv"
	n_synapses=10
	id_pre_layer=0
	id_post_layer=2
	output_dir="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/Networks/l1_to_l23"
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
	path_pre_h5="../Networks/nodes/layer_4_factor_0.1_nodes.h5"
	path_pre_csv="../Networks/nodes/layer_4_factor_0.1_node_types.csv"
	path_post_h5="../Networks/nodes/l5/layer_5_factor_0.1_nodes.h5"
	path_post_csv="../Networks/nodes/l5/layer_5_factor_0.1_node_types.csv"
	build_network_between_with_existing_nodes (path_pre_h5,path_pre_csv,path_post_h5,path_post_csv,df_connection_info_path,dict_types,n_synapses,output_dir)
	#build_network_between (dict_node_path,factor,df_connection_info_path,dict_types,n_synapses,id_pre_layer,id_post_layer,output_dir)

"""	
	
	
