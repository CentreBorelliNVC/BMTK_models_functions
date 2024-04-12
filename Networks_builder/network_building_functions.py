from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder
import plotly.express as px
import math
import numpy as np
import json
from position_building_V1 import get_population_loc_df,create_layer_population_loc,add_nodes_V1_in_nrrd
from edges_V1_functions import distance,distance_connection,distance_edges_within

def build_network (dict_node,factor,df_connection_info_path,dict_types,n_synapses,output_dir) : 
    net_layers,dataframes= add_nodes_V1_in_nrrd(dict_node_path,factor)
    v1=distance_edges_within(net_layers[2],df_connection_info_path,dict_types,n_synapses)
    v1.build()
    v1.save(output_dir)



if __name__ == '__main__':
	#node parameters
	dict_node_path="../Additional_data/dict_v1_nodes.json"
	factor=0.1
	
	#edge parameters
	df_connection_info=pd.read_csv("../Additional_data/l4_to_l4_connections_info.csv")
	dict_types=dict()
	dict_types["exc"]=["exc1","exc2","exc3","exc4","exc5","exc6","exc7"]
	dict_types["pvalb"]=["PV1","PV2"]
	dict_types["sst"]=["SST1","SST2","SST3"]
	dict_types["vip"]=["VIP1","VIP2","VIP3","VIP4"]
	dict_types["htr3a"]=["htr3a"]
	n_synapses=10
	
	#output parameter
	output_dir="../Networks/l4_to_l4"

	build_network (dict_node_path,factor,df_connection_info,dict_types,n_synapses,output_dir)
	
	
	
	
