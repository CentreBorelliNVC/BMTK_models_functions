#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 10:36:44 2024

@author: julienballbe
"""

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
from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
from bmtk.analyzer.spike_trains import plot_raster
from bmtk.analyzer.compartment import plot_traces
import json
import os
import plotly.express as px  
import math

import position_building_V1 as pos_build

Layer_pop_dict_L5 = {'Pyr5':.80,
                  "VIP" : .20,
                  "Total_N" : 2000,
                  'd_min' : 10,
                  'center_heigth' : -157.-572.-900.-(1411/2),
                  "height" : 1411.0,
                  "radius" : 200.0}

Layer_pop_dict_L6 = {'Pyr6':.80,
                  "SST" : .20,
                  "Total_N" : 2000,
                  'd_min' : 10,
                  'center_heigth' : -157.-572.-900.-1411.-(1973/2),
                  "height" : 1973.0,
                  "radius" : 200.0}

Full_layer_list = [Layer_pop_dict_L5, Layer_pop_dict_L6]



# test_dict_columnar={'granularity':1,
# 'radius':1000,
# 'd_min':10}

def create_network_geometry(geom_type, param_dict, layer_list):
    
    network_geom = CellLocations 
    
    
    if geom_type == "Column":
        if 'granularity' not in param_dict.keys or 'radius' not in param_dict.keys or 'd_min' not in param_dict.keys:
            raise ValueError("param_dict misses parameters")
            
        network_geometry = pos_build.create_column(param_dict, layer_list, do_plot = False)
        
        
        #L1-L2/3: 157 ± 16 μm, L2/3-L4: 575 ± 57 μm, L4-L5: 900 ± 50 μm, L5-L6: 1,411 ± 28 μm; L6-white matter (WM): 1,973 ± 44 μm
        
        print("test")
        
def get_nodes_csv_and_h5_table(network_name,network_folder):
    nodes_csv_table = pd.read_csv(f'{network_folder}/{network_name}_node_types.csv',sep=' ')
    
    
    current_file = h5py.File(f'{network_folder}/{network_name}_nodes.h5', 'r')
    
    
    current_file_node_group_id = pd.DataFrame(current_file['nodes'][network_name]['node_group_id'], columns=['node_group_id'])
    current_file_node_group_index = pd.DataFrame(current_file['nodes'][network_name]['node_group_index'], columns=['node_group_index'])
    current_file_node_id = pd.DataFrame(current_file['nodes'][network_name]['node_id'], columns=['node_id'])
    current_file_node_type_id = pd.DataFrame(current_file['nodes'][network_name]['node_type_id'], columns=['node_type_id'])
    test_df = pd.concat([current_file_node_group_id,current_file_node_group_index, current_file_node_id, current_file_node_type_id],axis=1)

    dynamics_params_list= current_file['nodes'][network_name]['0']['dynamics_params']
    dynamics_params_list_str = [x.decode('ascii') for x in dynamics_params_list]
    node_params_table = pd.DataFrame(dynamics_params_list_str,columns=['dynamics_params'])
    
    node_position_table = pd.DataFrame(current_file['nodes'][network_name]['0']['positions'],columns=['x','y','z'])
    
    node_params_position_table = pd.concat([test_df,node_params_table,node_position_table], axis=1)
    current_file.close()
    return nodes_csv_table, node_params_position_table
        


        
        
    
    


