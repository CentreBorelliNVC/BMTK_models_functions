#!/usr/bin/env python3
import json
import pandas as pd
import numpy as np
import sys
sys.path.append("../Networks_builder")
from network_building_functions import build_network

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

