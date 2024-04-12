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


	
	
	
	
