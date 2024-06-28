#!/usr/bin/env python3
from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder
import plotly.express as px
import math
import numpy as np
import json
import tqdm
from position_building_V1 import get_population_loc_df,create_layer_population_loc,add_nodes_V1_in_nrrd
import General_functions as gen_func

def distance (x_pre,x_post,z_pre,z_post) :
	d=np.sqrt((x_pre-x_post)**2+(z_pre-z_post)**2)
	return(d)

def delay_function (source,target,velocity) : #velocity in µm/ms 
	d=distance(source['positions'][0],target['positions'][0],source['positions'][2],target['positions'][2])
	delay=round(d/velocity,1) #the delay cannot be smaller than the dt (by default dt = 0.1ms)
	if delay < 0.1  :
		delay =0.1
	return(delay)

def velocity_extraction(velocity_path,layer_pre,layer_post) :
	df=pd.read_csv(velocity_path)
	velocity=list(df.loc[df["Unnamed: 0"]==layer_pre][layer_post])[0]
	return(velocity)
	
def distance_orientation_connection(source,target,amplitude,mu,function_type,n_synapses,slope) : #only works for now for nodes h5 with positions (instead of x,y,z) 
	
	#print(source['location'],source['pop_name'],source['node_id'],target["location"],target['pop_name'],target['node_id'])
	d=distance(source['positions'][0],target['positions'][0],source['positions'][2],target['positions'][2])
	if function_type == "gaussian" : 
		proba_d=amplitude*np.exp(-d**2/(2*mu**2))
	else : 
		if function_type=="exponential" : 
			proba_d=amplitude*np.exp(-d/mu)
		else : 
			proba_d=amplitude*np.exp(-d**2/(2*mu**2)) #billeh function
	delta_orientation=float(source['tuning_angle'])-float(target['tuning_angle'])
	delta_orientation_norm = abs(abs(abs(180.0-abs(delta_orientation))-90.0)-90.0)
	p_angle=slope*delta_orientation_norm + 1
	if source["ei"]=="e" and target["ei"]=="e" :
		p=(proba_d/100)*p_angle
	else : 
		p=proba_d/100
	p_connect=np.random.binomial(1,p)
	if p_connect == 1:
		if type(n_synapses)!=int  : 
			if source["ei"]=="e" : 
				nsyns=np.random.randint(n_synapses[0][0],n_synapses[0][1])
			else : 
				nsyns=np.random.randint(n_synapses[1][0],n_synapses[1][1])
			print(nsyns)
			return(nsyns)
		else : 
			nsyns=n_synapses
			print(nsyns)
			return(nsyns)
	else : 
		print(0)
		return(0)	

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
		if type(n_synapses)!=int  :  
			if source["ei"]=="e" : 
				nsyns=np.random.randint(n_synapses[0][0],n_synapses[0][1])
			else : 
				nsyns=np.random.randint(n_synapses[1][0],n_synapses[1][1])
			print(nsyns)
			return(nsyns)
		else : 
			nsyns=n_synapses
			print(nsyns)
			return(nsyns)
	else : 
		return(0)


def distance_edges_within (net,df_connection_info_path,dict_types,n_synapses) : #faire une fonction qui choisit le nomber de n ou le mettre dans le dict
	df_connection_info=pd.read_csv(df_connection_info_path)
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

def distance_edges_within_orientation (net,df_connection_info_path,dict_types,n_synapses,slope) : #faire une fonction qui choisit le nomber de n ou le mettre dans le dict
	df_connection_info=pd.read_csv(df_connection_info_path)
	for i in np.arange(df_connection_info.shape[0]) : 
		pre_type=df_connection_info.loc[i]["pre"]
		post_type=df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[pre_type] : 
			for post_subtype in dict_types[post_type] : 
				net.add_edges(
					source={'pop_name':pre_subtype},
					target={'pop_name':post_subtype},
					connection_rule=distance_orientation_connection,
					connection_params={'amplitude':df_connection_info.loc[i]["pmax"],'mu':df_connection_info.loc[i]["sigma"],'function_type':df_connection_info.loc[i]["rule"],'n_synapses':n_synapses,'slope':slope},
        				syn_weight=df_connection_info.loc[i]["weight"],
					delay= np.random.uniform(1,3),
					dynamics_params=df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
	return(net)

def distance_orientation_custom_delay_within (net,df_connection_info_path,dict_types,n_synapses,slope,path_velocity,layer) : 
	df_connection_info=pd.read_csv(df_connection_info_path)
	velocity=velocity_extraction(path_velocity,layer,layer)
	for i in np.arange(df_connection_info.shape[0]) : 
		pre_type=df_connection_info.loc[i]["pre"]
		post_type=df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[pre_type] : 
			for post_subtype in dict_types[post_type] : 
				conn=net.add_edges(
					source={'pop_name':pre_subtype},
					target={'pop_name':post_subtype},
					connection_rule=distance_orientation_connection,
					connection_params={'amplitude':df_connection_info.loc[i]["pmax"],'mu':df_connection_info.loc[i]["sigma"],'function_type':df_connection_info.loc[i]["rule"],'n_synapses':n_synapses,'slope':slope},
        				syn_weight=df_connection_info.loc[i]["weight"],
					dynamics_params=df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
				conn.add_properties('delay',
					rule=delay_function,
					rule_params={'velocity': velocity},
					dtypes=float)
	return(net)



	
def distance_edges_within_with_custom_delay (net,df_connection_info_path,dict_types,n_synapses,path_velocity,layer) : #faire une fonction qui choisit le nomber de n ou le mettre dans le dict
	velocity=velocity_extraction(path_velocity,layer,layer)
	df_connection_info=pd.read_csv(df_connection_info_path)
	for i in np.arange(df_connection_info.shape[0]) : 
		pre_type=df_connection_info.loc[i]["pre"]
		post_type=df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[pre_type] : 
			for post_subtype in dict_types[post_type] : 
				conn=net.add_edges(
					source={'pop_name':pre_subtype},
					target={'pop_name':post_subtype},
					connection_rule=distance_connection,
					connection_params={'amplitude':df_connection_info.loc[i]["pmax"],'mu':df_connection_info.loc[i]["sigma"],'function_type':df_connection_info.loc[i]["rule"],'n_synapses':n_synapses},
        				syn_weight=df_connection_info.loc[i]["weight"],
					dynamics_params=df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
				conn.add_properties('delay',
					rule=delay_function,
					rule_params={'velocity': velocity},
					dtypes=float)				
	return(net)

def distance_edges_between_with_custom_delay (net_pre,net_post,df_connection_info_path,dict_types,n_synapses,path_velocity) :
	df=pd.read_csv(df_connection_info_path)
	name_pre=net_pre.name
	split_pre=name_pre.split("_")
	layer_pre="l"+str(split_pre[1])
	name_post=net_post.name
	split_post=name_post.split("_")
	layer_post="l"+str(split_post[1])
	velocity=velocity_extraction(path_velocity,split_pre[1],split_post[1])
	sub_df_connection_info=df.loc[(df["pre_layer"]==layer_pre)&(df["post_layer"]==layer_post)].dropna()
	sub_df_connection_info=sub_df_connection_info.drop(["Unnamed: 0"],axis=1)
	for i in list(sub_df_connection_info.index) : 
		pre_type=sub_df_connection_info.loc[i]["pre"]
		post_type=sub_df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[layer_pre][pre_type] : 
			for post_subtype in dict_types[layer_post][post_type] : 
				conn=net_pre.add_edges(
					source=net_pre.nodes(pop_name=pre_subtype),
					target=net_post.nodes(pop_name=post_subtype),
					connection_rule=distance_connection,
					connection_params={'amplitude':sub_df_connection_info.loc[i]["pmax"],'mu':sub_df_connection_info.loc[i]["sigma"],'function_type':sub_df_connection_info.loc[i]["rule"],'n_synapses':n_synapses},
        				syn_weight=sub_df_connection_info.loc[i]["weight"],
					dynamics_params=sub_df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
				conn.add_properties('delay',
					rule=delay_function,
					rule_params={'velocity': velocity},
					dtypes=float)
	
	return(net_pre)	


def distance_orientation_custom_delay_between(net_pre,net_post,df_connection_info_path,dict_types,n_synapses,slope,path_velocity) : 
	df=pd.read_csv(df_connection_info_path)
	name_pre=net_pre.name
	split_pre=name_pre.split("_")
	layer_pre="l"+str(split_pre[1])
	name_post=net_post.name
	split_post=name_post.split("_")
	layer_post="l"+str(split_post[1])
	velocity=velocity_extraction(path_velocity,split_pre[1],split_post[1])
	sub_df_connection_info=df.loc[(df["pre_layer"]==layer_pre)&(df["post_layer"]==layer_post)].dropna()
	sub_df_connection_info=sub_df_connection_info.drop(["Unnamed: 0"],axis=1)
	for i in list(sub_df_connection_info.index) : 
		pre_type=sub_df_connection_info.loc[i]["pre"]
		post_type=sub_df_connection_info.loc[i]["post"]
		for pre_subtype in dict_types[layer_pre][pre_type] : 
			for post_subtype in dict_types[layer_post][post_type] : 
				conn=net_pre.add_edges(
					source=net_pre.nodes(pop_name=pre_subtype),
					target=net_post.nodes(pop_name=post_subtype),
					connection_rule=distance_orientation_connection,
					connection_params={'amplitude':sub_df_connection_info.loc[i]["pmax"],'mu':sub_df_connection_info.loc[i]["sigma"],'function_type':sub_df_connection_info.loc[i]["rule"],'n_synapses':n_synapses,'slope':slope},
        				syn_weight=sub_df_connection_info.loc[i]["weight"],
					dynamics_params=sub_df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
				conn.add_properties('delay',
					rule=delay_function,
					rule_params={'velocity': velocity},
					dtypes=float)
	return(net_pre)



def distance_edges_between_orientation (net_pre,net_post,df_connection_info_path,dict_types,n_synapses,slope) : #faire une fonction qui choisit le nomber de n ou le mettre dans le dict (arugment devrait être un interval)
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
					connection_rule=distance_orientation_connection,
					connection_params={'amplitude':sub_df_connection_info.loc[i]["pmax"],'mu':sub_df_connection_info.loc[i]["sigma"],'function_type':sub_df_connection_info.loc[i]["rule"],'n_synapses':n_synapses,'slope':slope},
        				syn_weight=sub_df_connection_info.loc[i]["weight"],
					delay= np.random.uniform(1,3),
					dynamics_params=sub_df_connection_info.loc[i]["synaptic_type"],
					model_template='static_synapse'
				)
	return(net_pre)

	
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
	
def get_selection_probability(src_type, lgn_models_subtypes_dictionary):
  current_model_subtypes = lgn_models_subtypes_dictionary[src_type[0:4]]['sub_types']
  current_model_probabilities = lgn_models_subtypes_dictionary[src_type[0:4]]['probabilities']
  lgn_model_idx = [i for i, model in enumerate(current_model_subtypes) if src_type == model][0]
  return current_model_probabilities[lgn_model_idx]


def select_lgn_sources(sources, target, lgn_mean, probability, poissonParameter, sON_ratio, centers_d_min,
                       centers_d_max, ON_OFF_w_min, ON_OFF_w_max, aspectRatio_min, aspectRatio_max, N_syn):
  
  #print(target.node_id,target["pop_name"])
  #for each target we have all sources taken into account
  source_ids = [s.node_id for s in sources]
  # Check if target supposed to get a connection and if not, then no need to keep calculating.
  if np.random.random() > probability: #probability = 1 for all by default in the lgn_connection dict
    return [None] * len(source_ids)

  subfields_centers_distance_L = centers_d_max - centers_d_min
  subfields_ON_OFF_width_L = ON_OFF_w_max - ON_OFF_w_min
  subfields_width_aspect_ratio_L = aspectRatio_max - aspectRatio_min
  x_position_lin_degrees = np.tan(0.07 * np.array(target['positions'][0]) * np.pi / 180.0) * 180.0 / np.pi
  #x_position_lin_degrees = np.tan(0.07 * np.array(target['x']) * np.pi / 180.0) * 180.0 / np.pi
  y_position_lin_degrees = np.tan(0.04 * np.array(target['positions'][2]) * np.pi / 180.0) * 180.0 / np.pi
  #y_position_lin_degrees = np.tan(0.04 * np.array(target['z']) * np.pi / 180.0) * 180.0 / np.pi
  vis_x = lgn_mean[0] + ((x_position_lin_degrees))  # - l4_mean[0]) / l4_dim[0]) * lgn_dim[0]
  vis_y = lgn_mean[1] + ((y_position_lin_degrees))  # - l4_mean[2]) / l4_dim[2]) * lgn_dim[1]

  ellipse_center_x0 = vis_x  # tar_cells[tar_gid]['vis_x']
  ellipse_center_y0 = vis_y  # tar_cells[tar_gid]['vis_y']

  tuning_angle = float(target['tuning_angle'])
  tuning_angle = None if np.isnan(tuning_angle) else tuning_angle
  # tuning_angle = None if math.isnan(target['tuning_angle']) else target['tuning_angle']
  if tuning_angle is None:
    ellipse_b0 = (ON_OFF_w_min + np.random.uniform(0.0,
                                                   1.0) * subfields_ON_OFF_width_L) / 2.0  # Divide by 2 to convert from width to radius.
    ellipse_b0 = 2.5 * ellipse_b0  # 1.5 * ellipse_b0
    ellipse_a0 = ellipse_b0  # ellipse_b0
    top_N_src_cells_subfield = 15  # 20
    ellipses_centers_halfdistance = 0.0
    tuning_angle_value = 0.0
  else:
    tuning_angle_value = float(tuning_angle)
    ellipses_centers_halfdistance = (centers_d_min + np.random.uniform(0.0, 1.0) * subfields_centers_distance_L) / 2.0
    ellipse_b0 = (ON_OFF_w_min + np.random.uniform(0.0,
                                                   1.0) * subfields_ON_OFF_width_L) / 2.0  # Divide by 2 to convert from width to radius.
    ellipse_a0 = ellipse_b0 * (aspectRatio_min + np.random.uniform(0.0, 1.0) * subfields_width_aspect_ratio_L)
    ellipse_phi = tuning_angle_value + 180.0 + 90.0  # Angle, in degrees, describing the rotation of the canonical ellipse away from the x-axis.
    ellipse_cos_mphi = np.cos(-np.radians(ellipse_phi))
    ellipse_sin_mphi = np.sin(-np.radians(ellipse_phi))
    top_N_src_cells_subfield = 8  # 10 #9

    if np.random.random() < sON_ratio:
      cell_sustained_unit = 'sON_'
    else:
      cell_sustained_unit = 'sOFF_'

  cell_TF = np.random.poisson(poissonParameter)
  while cell_TF <= 0:
    cell_TF = np.random.poisson(poissonParameter)

  sON_subunits = np.array([1., 2., 4., 8.])
  sON_sum = np.sum(abs(cell_TF - sON_subunits))
  p_sON = (1 - abs(cell_TF - sON_subunits) / sON_sum) / (len(sON_subunits) - 1)

  sOFF_subunits = np.array([1., 2., 4., 8., 15.])
  sOFF_sum = np.sum(abs(cell_TF - sOFF_subunits))
  p_sOFF = (1 - abs(cell_TF - sOFF_subunits) / sOFF_sum) / (len(sOFF_subunits) - 1)

  tOFF_subunits = np.array([4., 8., 15.])
  tOFF_sum = np.sum(abs(cell_TF - tOFF_subunits))
  p_tOFF = (1 - abs(cell_TF - tOFF_subunits) / tOFF_sum) / (len(tOFF_subunits) - 1)
  #print(target.node_id,p_sON,p_sOFF,p_tOFF) #p_sON has 4 values ; p_sOFF 5 values ; p_tOFF 3 values
  # to match previous algorithm reorganize source cells by type
  cell_type_dict = {}
  lgn_models = ["sON_TF1","sON_TF2","sON_TF4", "sON_TF8","sOFF_TF1","sOFF_TF2","sOFF_TF4","sOFF_TF8","sOFF_TF15","tOFF_TF4","tOFF_TF8","tOFF_TF15","sONsOFF_001","sONtOFF_001"]
  for lgn_model in lgn_models:
    cell_type_dict[lgn_model] = [(src.node_id, src) for src in sources if src['model_name'] == lgn_model]
  lgn_models_subtypes_dictionary = {
    'sON_': {'sub_types': ['sON_TF1', 'sON_TF2', 'sON_TF4', 'sON_TF8'], 'probabilities': p_sON},
    'sOFF': {'sub_types': ['sOFF_TF1', 'sOFF_TF2', 'sOFF_TF4', 'sOFF_TF8', 'sOFF_TF15'],
             'probabilities': p_sOFF
             },
    'tOFF': {'sub_types': ['tOFF_TF4', 'tOFF_TF8', 'tOFF_TF15'], 'probabilities': p_tOFF},
  }

  src_cells_selected = {}
  for src_type in cell_type_dict.keys():
    src_cells_selected[src_type] = []

    if tuning_angle is None:
      ellipse_center_x = ellipse_center_x0
      ellipse_center_y = ellipse_center_y0
      ellipse_a = ellipse_a0
      ellipse_b = ellipse_b0
    else:
      if ('tOFF_' in src_type[0:5]):
        ellipse_center_x = ellipse_center_x0 + ellipses_centers_halfdistance * ellipse_sin_mphi
        ellipse_center_y = ellipse_center_y0 + ellipses_centers_halfdistance * ellipse_cos_mphi
        ellipse_a = ellipse_a0
        ellipse_b = ellipse_b0
      elif ('sON_' in src_type[0:5] or 'sOFF_' in src_type[0:5]):
        ellipse_center_x = ellipse_center_x0 - ellipses_centers_halfdistance * ellipse_sin_mphi
        ellipse_center_y = ellipse_center_y0 - ellipses_centers_halfdistance * ellipse_cos_mphi
        ellipse_a = ellipse_a0
        ellipse_b = ellipse_b0
      else:
        # Make this a simple circle.
        ellipse_center_x = ellipse_center_x0
        ellipse_center_y = ellipse_center_y0
        # Make the region from which source cells are selected a bit smaller for the transient_ON_OFF
        # cells, since each source cell in this case produces both ON and OFF responses.
        ellipse_b = ellipses_centers_halfdistance / 2.0  # 0.01
        ellipse_a = ellipse_b0  # 0.01 #ellipse_b0

    # Find those source cells of the appropriate type that have their visual space coordinates within the
    # ellipse.
    for src_id, src_dict in cell_type_dict[src_type]:
      x, y = (src_dict['x'], src_dict['y'])
      x = x - ellipse_center_x
      y = y - ellipse_center_y

      x_new = x
      y_new = y
      if tuning_angle is not None:
        x_new = x * ellipse_cos_mphi - y * ellipse_sin_mphi
        y_new = x * ellipse_sin_mphi + y * ellipse_cos_mphi

      if ((x_new / ellipse_a) ** 2 + (y_new / ellipse_b) ** 2) <= 1.0:
        #print("yes 1 : ((x_new / ellipse_a) ** 2 + (y_new / ellipse_b) ** 2) <= 1.0")
        if tuning_angle is not None:
          #print("tuning not none")
          if src_type == 'sONsOFF_001' or src_type == 'sONtOFF_001':
            #print("if yes 1 & source type sONsOFF or sONtOFF")
            src_tuning_angle = float(src_dict['tuning_angle'])
            delta_tuning = abs(abs(abs(180.0 - abs(tuning_angle_value - src_tuning_angle) % 360.0) - 90.0) - 90.0)
            if delta_tuning < 15.0:
              #print("if yes & source type sONsOFF or sONtOFF & delta_tuning < 15.0")
              src_cells_selected[src_type].append(src_id)

          # elif src_type in ['sONtOFF_001']:
          #     src_cells_selected[src_type].append(src_id)

          elif cell_sustained_unit in src_type[:5]:
            selection_probability = get_selection_probability(src_type, lgn_models_subtypes_dictionary)
            #print("if not sONsOFF/sONtOFF but is sON/sOFF, selected proba : ",selection_probability)
            if np.random.random() < selection_probability:
              #print("if is sON/sOFF, selected proba > random with source type :",src_type)
              src_cells_selected[src_type].append(src_id)

          elif 'tOFF_' in src_type[:5]:
            selection_probability = get_selection_probability(src_type, lgn_models_subtypes_dictionary)
            #print("if not sONsOFF/sONtOFF but is tOFF, selected proba is : ",selection_probability)
            if np.random.random() < selection_probability:
              #print("if not sONsOFF/sONtOFF but is tOFF, selected proba is > random")
              src_cells_selected[src_type].append(src_id)


        else:
          if (src_type == 'sONsOFF_001' or src_type == 'sONtOFF_001'):
            src_cells_selected[src_type].append(src_id)
            #print("if tuning angle is nan : ",tuning_angle," and is sONsOFF/sONtOFF")
          else:
            #print("lol")
            selection_probability = get_selection_probability(src_type, lgn_models_subtypes_dictionary)
            #print("if tuning angle is nan : ",tuning_angle," and is sON/sOFF/tOFF ; selected proba : ",selection_probability)
            if np.random.random() < selection_probability:
              #print("if tuning angle is nan : ",tuning_angle,"and is sON/sOFF/tOFF and selected proba > random")
              src_cells_selected[src_type].append(src_id)

        #print(src_cells_selected)
  select_cell_ids = [id for _, selected in src_cells_selected.items() for id in selected]
  #print(target.node_id,select_cell_ids)
  #a=[print(N_syn) if id in select_cell_ids else None for id in source_ids]
  return [N_syn if id in select_cell_ids else None for id in source_ids]

def lgn_to_v1_layer (net_pre,net_post,lgn_to_l4_dict,layer_types,field_size) :
	for pop_name in layer_types :
		#print(lgn_to_l4_dict.keys())
		conn_props = lgn_to_l4_dict[pop_name]['connection_params']
		conn_props['lgn_mean'] = (field_size[0]/2.0, field_size[1]/2.0)
		edge_props = lgn_to_l4_dict[pop_name]['edge_types_params']
		net_post.add_edges(
			source=net_pre.nodes(),
			target=net_post.nodes(pop_name=pop_name),
			connection_rule=select_lgn_sources,
			connection_params=conn_props,
			iterator='all_to_one',
			**edge_props
		)
	return(net_post)	

def add_edges_test(pre_net_name, pre_net_folder, post_net_name, post_net_folder, connection_table, n_synapses, slope, velocity_path, custom_delay):
    
    pre_net_csv, pre_net_h5 = gen_func.get_nodes_csv_and_h5_table(pre_net_name, pre_net_folder)
    post_net_csv, post_net_h5 = gen_func.get_nodes_csv_and_h5_table(post_net_name, post_net_folder)
    
    pre_net_pop_list = list(pre_net_csv.loc[:,'pop_name'])
    post_net_pop_list = list(post_net_csv.loc[:,'pop_name'])
    
    pre_net_obj = NetworkBuilder(name=pre_net_name)
    pre_net_obj.import_nodes(nodes_file_name=f"{pre_net_folder}/{pre_net_name}_nodes.h5", 
                             node_types_file_name=f"{pre_net_folder}/{pre_net_name}_node_types.csv")
    
    post_net_obj = NetworkBuilder(name=post_net_name)
    post_net_obj.import_nodes(nodes_file_name=f"{post_net_folder}/{post_net_name}_nodes.h5", 
                             node_types_file_name=f"{post_net_folder}/{post_net_name}_node_types.csv")

    established_connections = pd.DataFrame(columns = connection_table.columns)
    for pre_net_pop in tqdm.tqdm(pre_net_pop_list):

        #pre_layer, pre_pop = pre_net_pop.split("_") # To be modified when layer will be written L_2/3, L_4..
        pre_L, pre_layer_nb, pre_pop = pre_net_pop.split("_")
        pre_layer = f'{pre_L}_{str(pre_layer_nb)}'
        
        source_ei= np.array(pre_net_csv.loc[pre_net_csv['pop_name']==pre_net_pop,"ei"])[0]
        
        for post_net_pop in post_net_pop_list:
            post_L, post_layer_nb, post_pop = post_net_pop.split("_")
            post_layer = f'{post_L}_{str(post_layer_nb)}'
            
                
            sub_connection_table = connection_table.loc[(connection_table["Pre_Layer"]==pre_layer)&
                                                        (connection_table["Pre_pop"]==pre_pop)&
                                                        (connection_table["Post_Layer"]==post_layer)&
                                                        (connection_table["Post_pop"]==post_pop),:]

            if sub_connection_table.shape[0]==0:
                #print(f'No connection from {pre_net_pop} to {post_net_pop}')
                continue
            
            
            else:
                #Check if this connection has already been implemented
                sub_established_connections = established_connections.loc[(established_connections["Pre_Layer"]==pre_layer)&
                                                            (established_connections["Pre_pop"]==pre_pop)&
                                                            (established_connections["Post_Layer"]==post_layer)&
                                                            (established_connections["Post_pop"]==post_pop),:]
                if sub_established_connections.shape[0]!=0:
                    # if source to target connection has already been implemented do not re-do it
                    #print(f'Connection from {pre_net_pop} to {post_net_pop} already done')
                    
                    continue
                else:
                    pre_adapted_layer_number = pre_layer.split('_')[1]
                    post_adapted_layer_number = post_layer.split('_')[1]
                    
                    velocity=velocity_extraction(velocity_path, pre_adapted_layer_number, post_adapted_layer_number)
                    
                    established_connections = pd.concat([established_connections,sub_connection_table], ignore_index = True)
                    sub_connection_table = sub_connection_table.reset_index(drop=True)
                    
                    conection_rule = sub_connection_table.loc[0,"rule"]
                    connection_pmax = sub_connection_table.loc[0,"pmax"]
                    connection_sigma = sub_connection_table.loc[0,'sigma']
                    connection_weitgh = sub_connection_table.loc[0,"weight"]
                    connection_type = sub_connection_table.loc[0,"synaptic_type"]
                    
                    conn=pre_net_obj.add_edges(
    					source=pre_net_obj.nodes(pop_name=pre_net_pop),
    					target=post_net_obj.nodes(pop_name=post_net_pop),
    					connection_rule=distance_connection_test,
    					connection_params={'amplitude':connection_pmax,
                            'mu':connection_sigma,
                            'function_type':conection_rule,
                            'n_synapses':n_synapses,
                            "slope":None}, #if need orientation, enter not None slope
            			syn_weight=connection_weitgh,
    					dynamics_params=connection_type,
    					model_template='static_synapse'
    				)
                    conn.add_properties("delay",
                                               rule = delay_function_test,
                                               rule_params = {'velocity' : velocity,
                                                              'custom_delay' : custom_delay},# if need custom delay, set custom delay to True
                                               dtypes = float)
                    
             
    pre_net_obj.build()
    pre_net_obj.save_edges(output_dir = pre_net_folder)
    
    return established_connections

def delay_function_test(source, target, velocity, custom_delay):
    if custom_delay == False:
        delay = np.random.uniform(1,3)
    elif custom_delay == True:
        d=distance(source["positions"][0], target['positions'][0], source['positions'][2], target['positions'][2])
        delay = round(d/velocity, 1)
        if delay <0.1:
            delay=0.1
    return delay


def distance_connection_test(source, target, amplitude, mu, function_type, n_synapses, slope=None, orientation = False):
    d=distance(source['positions'][0], target['positions'][0], source["positions"][2], target["positions"][2])
    
    if function_type == 'gaussian':
        proba = amplitude*np.exp(-d**2/(2*mu**2))
    elif function_type == "exponential":
        proba = amplitude*np.exp(-d/mu)
    elif function_type == "billeh":
        proba = amplitude*np.exp(-d**2/(2*mu**2))
        
    if slope is not None: #In case we consider orientation 
        delta_orientation = float(source['tuning_angle'])-float(target['tuning_angle'])
        delta_orientation_norm = abs(abs(abs(180.0 - abs(delta_orientation))-90.0)-90.0)
        p_angle = slope*delta_orientation_norm + 1
        
        if source["ei"]=="e" and target["ei"]=='e':
            p = (proba/100)*p_angle
        else:
            p = proba/100
        p_connect = np.random.binomial(1/p)

    else: #in case we don't consider orientation
        p_connect = np.random.binomial(1, proba/100)

        
    if p_connect == 1:
        if type(n_synapses) != int :
            if source["ei"]=="e":
                nsyns = np.random.randint(n_synapses[0][0], n_synapses[0][1])
            else:
                nsyns = np.random.randint(n_synapses[1][0], n_synapses[1][1])
                
            return nsyns
        else:
            nsyns = n_synapses
            return nsyns
    else:
        return(0)

"""
if __name__ == '__main__':
	dict_path="../Additional_data/dict_v1_nodes.json"
	#net_layers,dataframes=add_nodes_V1_in_nrrd (dict_path,1) 
	df_connection_info_path="../Additional_data/l23_to_l23_connections_info.csv"
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
	path_velocity="/home/margaux/miniconda3/envs/ENV_NEST2/BMTK_models_functions/Additional_data/velocity.csv"
	path_pre_h5="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/Networks/nodes/layer_4_factor_0.1_nodes.h5"
	path_post_h5="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/Networks/nodes/layer_4_factor_0.1_nodes.h5"
	path_pre_csv="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/Networks/nodes/layer_4_factor_0.1_node_types.csv"
	path_post_csv="/home/margaux/miniconda3/envs/ENV_NEST2/stockage_github/Networks/nodes/layer_4_factor_0.1_node_types.csv"
	path_pre_split=path_pre_h5.split("/")
	path_pre_split_bis=path_pre_split[-1].split("_nodes")
	net_pre=NetworkBuilder(path_pre_split_bis[0])
	path_post_split=path_post_h5.split("/")
	path_post_split_bis=path_post_split[-1].split("_nodes")
	net_post=NetworkBuilder(path_post_split_bis[0])
	net_pre.import_nodes(nodes_file_name=path_pre_h5,node_types_file_name=path_pre_csv)
	net_post.import_nodes(nodes_file_name=path_post_h5,node_types_file_name=path_post_csv)
	#net=distance_edges_within(net_layers[2],df_connection_info,dict_types,n_synapses)
	#net=distance_edges_between(net_layers[0],net_layers[1],df_connection_info_path,dict_types,n_synapses) #changer argument
	net=distance_edges_between_with_custom_delay(net_pre,net_post,df_connection_info_path,dict_types,n_synapses,path_velocity)
	#net.build() #all connections within l4
	"""


	






