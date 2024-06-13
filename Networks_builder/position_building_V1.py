#!/usr/bin/env python3
from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder
import plotly.express as px
import math
import numpy as np
import json
import random
from random import sample
import copy
import General_functions as gen_func
import h5py
import warnings
import scipy
from scipy.stats import skewnorm, lognorm, norm
from scipy.special import erfc, erf
import plotnine as p9

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


#%%

def create_column(param_dict, layer_list, do_plot = False):
    '''
    Create a network cell location object. Given the parameters entered in the params dict, and the layer wanted, the function automatically computes the depth of the layer, the density, and the width
    
    Parameters
    ----------
    param_dict : Dictionnary
        "d_min" (float) representing the minimum distance in um between the nodes
        "radius" (float) representing the raidus in um of the column
        "granularity" (float) representing the proportion of biological cell density to keep in the column. Must be between 0 and 1.
    layer_list : List
        Contains the layers wanted in the columns. The layers must be designated as : L_1, L_2/3....
    do_plot : Bool, optional
        Open a web tab with a plotly representation of the cell positions. The default is False.

    Returns
    -------
    network_cell_location : CellLocations object
        DESCRIPTION.

    '''
    
    Layer_width_dict = {"L_1" : 157.0,# in um
                        "L_2/3" : 575.0, 
                        "L_4" : 900.0,
                        "L_5" : 1411.0, 
                        "L_6" : 1973}
    
    Layer_center_dict = {"L_1" : -78.5, #in um
                         "L_2/3" : -444.5, 
                         "L_4" : -1182.0,
                         "L_5" : -2334.5,
                         "L_6" : -4026.5}
    
    Layer_density_dict = {"L_1" : 43_199.7, #from BBP, in /mm3
                         "L_2/3" : 108_787.9, 
                         "L_4" : 130_640.0,
                         "L_5" : 105_268.4,
                         "L_6" : 98_393.2}
    
    
    
    network_cell_location = CellLocations('ctx')
        
    for current_layer in layer_list:
        population_name = current_layer
        #compute desired number of cell in the layer, according to the radius of the cynlinder, and the granularity
        granularity = param_dict['granularity']

        radius = param_dict['radius']
        d_min = param_dict['d_min']
        
        
        cylinder_volume = math.pi * radius**2 * Layer_width_dict[population_name] # in um3
        N = round(cylinder_volume * 
                  (Layer_density_dict[population_name]*1e-9) * # convert density in um
                  granularity)
        
        
        
        network_cell_location.dmin = d_min
        
        
        network_cell_location.add_positions_columnar(pop_names=population_name,
                                   N=N, 
                                   center=[0.0, Layer_center_dict[population_name] , 0.0], 
                                   height = Layer_width_dict[population_name],
                                   min_radius = 0.0, 
                                   max_radius = param_dict['radius'], 
                                   method='prog', 
                                   verbose=False)
        
    position_df = get_population_loc_df(network_cell_location)
    if do_plot:
        plot_population(position_df)

    
    return network_cell_location

def plot_population(position_df):
    
    fig = px.scatter_3d(position_df, x='x', y='z', z='y',color='Cell_type')
    fig.update_traces(marker_size=4)
    fig.update_scenes(zaxis_autorange="reversed")
    fig.show('browser')
    
    return position_df

def plot_network(network_name,network_folder):
    
    network_csv_table, network_h5_table = gen_func.get_nodes_csv_and_h5_table(network_name,network_folder)
    network_h5_table = network_h5_table.astype({'node_type_id':'int'})
    Full_network_table = pd.merge(network_csv_table,network_h5_table, on ='node_type_id')
    
    fig = px.scatter_3d(Full_network_table, x='x', y='z', z='y',color='pop_name')
    fig.update_traces(marker_size=4)
    fig.update_scenes(zaxis_autorange="reversed")
    fig.show('browser')
    
    
    
def get_population_loc_df(cell_location_obj):
    # /!\ added a condition for the tuning angle (if tuning_angle = True --> add tuning angle column)
    for pop_name in cell_location_obj._all_pop_names : 
    	if hasattr(getattr(cell_location_obj,pop_name),"tuning_angles") == True :
    		column_names=["x","y","z","Tuning_angle","Cell_type"]
    	else : 
    		column_names=["x","y","z","Cell_type"]
    position_df=pd.DataFrame(columns=column_names)
    for pop_name,pop_loc in zip(cell_location_obj._all_pop_names,cell_location_obj._all_positions):
        
        pop_loc_df=pd.DataFrame(pop_loc)
        pop_loc_df.columns = ['x','y','z']
        if "Tuning_angle" in column_names : 
        	pop_loc_df.loc[:,"Tuning_angle"] = list(getattr(cell_location_obj,pop_name).tuning_angles)
        pop_loc_df.loc[:,'Cell_type'] = pop_name
        position_df = pd.concat([position_df, pop_loc_df],ignore_index=True)
    return position_df
        
        

def create_layer_population_loc(Layer_pop_dict, ctx=None, do_plot=False):
    geom_construction_params = ['d_min','center_heigth','height','radius','Total_N']
    Sub_Layer_pop_dict = {key: value for key, value in Layer_pop_dict.items() if key not in geom_construction_params}
    
    pop_name_list = list(Sub_Layer_pop_dict.keys())
    pop_prop_list =list(Sub_Layer_pop_dict.values())
    
    if ctx == None:
        ctx = CellLocations('ctx')
        
    
    ctx.dmin = Layer_pop_dict['d_min']
    
    ctx.add_positions_columnar(pop_name_list ,
                               partitions=pop_prop_list, 
                               N=Layer_pop_dict['Total_N'], 
                               center=[0.0, Layer_pop_dict['center_heigth'], 0.0], 
                               height = Layer_pop_dict['height'],
                               min_radius = 0.0, 
                               max_radius = Layer_pop_dict['radius'], 
                               method='prog', 
                               verbose=False)
    total_pop = 0
    
    # compute density for each population and total density
    for pop_name,pop_loc in zip(ctx._all_pop_names,ctx._all_positions):
        pop_size = pd.DataFrame(pop_loc).shape[0]
        total_pop+=pop_size
        cylinder_volume = math.pi * Layer_pop_dict['radius']**2 * Layer_pop_dict['height']
        population_density = pop_size/cylinder_volume
        print(f'Density of population {pop_name} = {population_density} cells/mm3')
        
    total_density = total_pop / cylinder_volume
    print(f'Density of total population = {total_density} cells/mm3')
    
    
    position_df = get_population_loc_df(ctx)
    if do_plot:
        plot_population(position_df)

    return ctx, position_df

def select_models(model_feature_table, N ,method = 'Uniform', fit_table = None, feature = None, population=None, layer=None, do_plot = False):
    '''
    Generate a list of N models, according to a selection method. 
    Uniform will take all available models the same amount of time
    Distributed selects models so that the models' feature distirbution follows the corresponding distribution indicated in fit_table
    Average selects one model whic is the closest to the feature average indicated in fit_table

    Parameters
    ----------
    model_feature_table : pd.DataFrame
        Table gathering for all available models the different properties.
    N : Int
        Number of models to select.
    method : str, optional
        Specify the method to delect the models. Must be Uniform, Distributed or Average. The default is 'Uniform'.
    fit_table : pd.DataFrame, optional
        If method not Uniform, use fit table to have the feature distribution fit parameters per population per layer. The default is None.
    feature : str, optional
        If method not Uniform, feature of interest for which the distribution must be respected . The default is None.
    population : str, optional
        If method not Uniform, population of interest for which the distribution must be respected . The default is None.
    layer : str, optional
        If method not Uniform, layer of interest for which the distribution must be respected . The default is None.
    do_plot : Boolean, optional
        if method is 'Distributed' plot the selected models' features distribution against the theorical feature distribution. The default is False.

    Returns
    -------
    selected_cell_type_model : List
        List of models file selected.
    selected_cell_type_model_df : pd.DataFrame
        DataFrame of selected models file with models' features .

    '''
    
    selected_cell_type_model = []
    selected_cell_type_model_df = pd.DataFrame(columns = model_feature_table.columns)

    if method == 'Uniform':
        
        # select desired number of model file
        # All models are selected enough time too have n models 
        model_list = list(model_feature_table.loc[:,'file_name'])
        while len(selected_cell_type_model) < N: 
            
            needed_model_id = N-len(selected_cell_type_model)
            
            # here each model has an equal chance of being selected --> uniform distribution 
            current_selected_models = sample(model_list, min(N,needed_model_id,len(model_list) ))
            selected_cell_type_model = selected_cell_type_model + current_selected_models
    
    elif method == 'Distributed':
        # Given a feature distribution, attribute to each model for which the feature has been computed, a probability of being selected
        # the probability depends on the the observed feature distribution 
        #model_feature_table = model_feature_table.dropna([feature], inplace=True)
        population_fit_line = fit_table.loc[(fit_table['Population']==population)&(fit_table['Layer_second']==layer)&(fit_table['Feature']==feature),:]

        population_fit_line=population_fit_line.reset_index(drop=True)
        fit = population_fit_line.loc[0,'Fit']
        A = population_fit_line.loc[0,'A']
        mu = population_fit_line.loc[0,'mu']
        sigma = population_fit_line.loc[0,'sigma']
        gamma = population_fit_line.loc[0,'gamma']
        bin_width =population_fit_line.loc[0,'Bin_width']
        
        distrib_value_list = get_list_from_distrib(fit, N, mu, sigma, gamma)

        feature_value = np.array(model_feature_table.loc[:,feature])
        for val in distrib_value_list:
            closest_model_val = min(feature_value, key=lambda x: abs(x - val))
            closest_model_row = model_feature_table.loc[model_feature_table[feature] == closest_model_val,:]
            if closest_model_row.shape[0]!=1: # in case two models are equally close to desired number
                closest_model_row = closest_model_row.sample(n=1)
                
            
            selected_cell_type_model_df = pd.concat([selected_cell_type_model_df, closest_model_row], ignore_index=True)
        

        selected_cell_type_model = list(selected_cell_type_model_df.loc[:,'file_name'])
        
        if do_plot:
            min_x = np.nanmin(selected_cell_type_model_df.loc[:,feature])
            max_x = np.nanmax(selected_cell_type_model_df.loc[:,feature])
            amp_x = max_x - min_x
            test_x_data=np.arange(min_x-.1*amp_x, max_x+.1*amp_x, .01*amp_x)
            if fit == 'Skewed Gaussian':
                theorical_distrib = skewedgaussian(test_x_data, A, gamma, mu, sigma)
            elif fit == "Gaussian":
                theorical_distrib = gaussian(test_x_data, A, mu, sigma)
            elif fit == "LogNormal":
                theorical_distrib = Log_normal(test_x_data, A, mu, sigma)
                
            
            theorical_ditrib_df = pd.DataFrame({'Feature':test_x_data, 'Count':theorical_distrib})
            
            bin_df = get_data_bin_df(np.array(selected_cell_type_model_df.loc[:,feature]), bin_width)

            multiplicativ_factor = np.nanmax(bin_df.loc[:,'Count'])/np.nanmax(theorical_ditrib_df.loc[:,'Count'])
            theorical_ditrib_df.loc[:,'Count'] *= multiplicativ_factor
            
            my_plot = p9.ggplot(selected_cell_type_model_df, p9.aes(x=feature))+p9.geom_histogram(binwidth = bin_width)
            my_plot += p9.geom_line(theorical_ditrib_df, p9.aes(x="Feature", y="Count"))
            my_plot += p9.ggtitle(f'{population} neurons in layer {layer} , {feature} distribution')
            my_plot.show()
            
    elif method == "Average":
        #For each layer-population, take the model whose feature is the closest to the distribution mean
        
        population_fit_line = fit_table.loc[(fit_table['Population']==population)&(fit_table['Layer_second']==layer)&(fit_table['Feature']==feature),:]
        population_fit_line=population_fit_line.reset_index(drop=True)

        mu = population_fit_line.loc[0,'mu']
        
        feature_value = np.array(model_feature_table.loc[:,feature])
        closest_model_val = min(feature_value, key=lambda x: abs(x - mu))
        closest_model_row = model_feature_table.loc[model_feature_table[feature] == closest_model_val,:]

        for elt in range(N):
            
            selected_cell_type_model_df = pd.concat([selected_cell_type_model_df, closest_model_row], ignore_index=True)
        
        selected_cell_type_model = list(selected_cell_type_model_df.loc[:,'file_name'])
        
    return selected_cell_type_model, selected_cell_type_model_df
        
        

def get_list_from_distrib(model,N,  mu, sigma, gamma):
    """
    Generate for a given function a list of N random numbers

    Parameters
    ----------
    model : str
        Gaussian, LogNormal or Skeweed Gaussian.
    N : int
        number of float to generate.
    mu : float

    sigma : float

    gamma : float


    Returns
    -------
    list_values : list
        list of float following the desired distribution.

    """
    if model =='Skewed Gaussian':
        list_values = skewnorm.rvs( a=gamma, loc=mu, scale=sigma, size=N)

    elif model == 'LogNormal':
        list_values = lognorm.rvs(size= N, s=sigma, scale=np.exp(mu))
        
    elif model == 'Gaussian':
        list_values = norm.rvs(size=N, loc=mu, scale=sigma)
        
    return list_values


def skewedgaussian(x, A, gamma, mu, sigma):
    y = (A / (sigma * np.sqrt(2 * np.pi))) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) * (
                1 + erf((gamma * (x - mu)) / (sigma * np.sqrt(2))))
    return (y)


def get_data_bin_df(value_array, bin_width_or_nb, bin_def = "Width"):
    
    if bin_def == 'Number':
        bin_edges = np.linspace(np.nanmin(value_array), np.nanmax(value_array), bin_width_or_nb + 1) # includes right edge, so adding one to desired bin number
        bin_width = bin_edges[1] - bin_edges[0]
        bin_edges = np.linspace((np.nanmin(value_array)-bin_width), (np.nanmax(value_array)+bin_width),bin_width_or_nb)
        
    else:

        bin_nb = int(((np.nanmax(value_array)+bin_width_or_nb) - (np.nanmin(value_array)-bin_width_or_nb))//bin_width_or_nb)

        bin_edges = np.linspace((np.nanmin(value_array)-bin_width_or_nb), (np.nanmax(value_array)+bin_width_or_nb),bin_nb) # includes right edge, so adding one to desired bin number
        
        
    current_output,current_bin_edges,current_bins_number = scipy.stats.binned_statistic(value_array,
                                    value_array,
                                    statistic='count',
                                    bins=bin_edges)

    bin_center_array=current_bin_edges[1:]-(current_bin_edges[1:]-current_bin_edges[:-1])/2
    bin_count_df=pd.DataFrame({"Feature":bin_center_array,"Count":current_output})
    
    return bin_count_df

def gaussian(x, A, mu, sigma):
    y = (A / (sigma * np.sqrt(2 * np.pi))) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
    return (y)
            
def Log_normal(x, A, mu, sigma):
    y = (A / (sigma * np.sqrt(2 * np.pi))) * (np.exp(-((np.log(x) - mu) ** 2) / (2 * (sigma ** 2))) / x)
   
    return (y)

def function_create_nodes_dict(model_features_table, layer_prop_dict, position_obj, Model_population_table, feature_distribution_table):
    '''
    Given a dictionnary describing for each layer, the populations their proportion, and how to select the models for eahc of the population
    

    Parameters
    ----------
    model_features_table : pd.DataFrame
        Table gathering for all available models the different properties.
    layer_prop_dict : Dict
        Dictionnary describing for each layer, the populations their proportion, and how to select the models for eahc of the population .
        layer_prop_dict = {'L_4':{'PV':{'Proportion':.8,
                                       "Model_selection":"Average",
                                       "Feature":'Time_constant_ms'}, 
                                 "Sst":{'Proportion':.2,
                                        "Model_selection":"Uniform"}}}
    position_obj 
        BMTK object containing x,y,z coordintes for each layer. Generated by create_column function.
    Model_population_table : pd.DataFrame
        Table containing for each cell the different categroies they fall into. Notably Exc_Inh, Exc_Inh_FT and Cell_type which must be present in the table .
    feature_distribution_table : pd.DataFrame
        Table gathering feature distribution fit parameters per population per layer.

    Raises
    ------
    ValueError
        If the layer description dictionnary misses some information, raises Value Error.

    Returns
    -------
    node_dict : Dictionnary
        Dict that will be used to generate the nodes of teh network, with appropriate population proportion, coordinates and models.

    '''

    node_dict={}
    N_tot = 0
    
    #Make sure the model_id are str
    Model_population_table  = Model_population_table.dropna(subset=['Model_id'])
    Model_population_table = Model_population_table.astype({"Model_id":"int"})
    Model_population_table = Model_population_table.astype({"Model_id":"str"})
    
    model_features_table  = model_features_table.dropna(subset=['Model_id'])
    model_features_table = model_features_table.astype({"Model_id":"int"})
    model_features_table = model_features_table.astype({"Model_id":"str"})
    
    
    
    #Check that layer_prop_dict has all needed information
    for Layer, Layer_dict  in layer_prop_dict.items():
        layer_prop=0
        for pop, pop_dict in Layer_dict.items():
            layer_prop+=pop_dict['Proportion']
            if pop_dict['Model_selection'] != "Uniform" and "Feature" not in pop_dict.keys():
                raise ValueError(f'For Layer {Layer}, population {pop}, missing feature as Model selection not Uniform')
                
            if pop not in Model_population_table.loc[:,"Exc_Inh"].unique() and pop not in Model_population_table.loc[:,"Exc_Inh_FT"].unique() and pop not in Model_population_table.loc[:,"Cell_type"].unique():
                raise ValueError(f'For Layer {Layer}, population {pop} not in model population table')
        if np.round(layer_prop,3) != 1.:
            raise ValueError(f'For Layer {Layer}, total proportion must be 1 (here = {np.round(layer_prop,3)})')
        
    # for each layer, get model list to respect specified cell type proportion 
    population_position_df = get_population_loc_df(position_obj)
    
   
    for layer, layer_pop_dict in layer_prop_dict.items(): 
        position_obj_modified = copy.deepcopy(position_obj)
        layer_positions_list = position_obj_modified._all_pop_names

        # get the position array that correspond to the layer of interest, so that we can partition it in the cell type loop
        layer_index = layer_positions_list.index(layer)
        layer_position_array = position_obj_modified._all_positions[layer_index]

        
        # split coordinates so that we can pass them direclty to add nodes cell-type-wise
        Layer_positions_df = population_position_df.loc[population_position_df.loc[:,'Cell_type'] == layer,:]
        layer_N = Layer_positions_df.shape[0] # get number of position available for the current layer
        Layer_positions_df = Layer_positions_df.reset_index(drop=True)
        Layer_positions_index = Layer_positions_df.index.values
        
        # randomly partition position_df indexes
        random.shuffle(Layer_positions_index)
        
        
        # for each cell type get model list to respect specified cell type proportion 
        last_index = 0
        
        layer_node_dict = {'Layer' : layer,
                           'N' : layer_N}

        
        layer_node_count={}
        total_layer = 0
        # determine number of nodes required for this cell type
        for population, population_in_layer_dict in layer_pop_dict.items():
            cell_type_in_layer_proportion = population_in_layer_dict['Proportion']
            cell_type_nb = int(round(cell_type_in_layer_proportion * layer_N))
            layer_node_count[population]=cell_type_nb
            total_layer += int(round(cell_type_in_layer_proportion * layer_N))
        #Check if the addition of the populations' sizes add up to the layer size, if not, adjust on the most represented population
        if total_layer != layer_N:
            diff = layer_N - total_layer
            
            max_key = next(iter(layer_node_count))
            for key in layer_node_count:
                # If the value of the current key is greater than the value of max_key, update max_key
                if layer_node_count[key] > layer_node_count[max_key]:
                    max_key = key
            
            #Add to the most represented population the difference, in order to minimize the relative proportion added
            layer_node_count[max_key] += diff
                
        
        for population, population_in_layer_dict in layer_pop_dict.items():
            # create a copy of original position obj so we can modify it
            cell_type_position_obj = copy.deepcopy(position_obj)
            
            cell_type_pop_name = f'{layer}_{population}'
            
            # is the cell type exc or inh
            if 'Exc'.casefold() in population.casefold():
                ei='e'
            else:
                ei='i'
            
            # determine number of nodes required for this cell type
            cell_type_in_layer_proportion = population_in_layer_dict['Proportion']
            
            desired_cell_Type_count = layer_node_count[population]
            N_tot += desired_cell_Type_count
            
            # get random list of coordinated, one coord for each node, and select the corresponding list
            layer_cell_type_position_index =  Layer_positions_index[last_index:(last_index+desired_cell_Type_count)]
            layer_cell_type_position_array = layer_position_array[layer_cell_type_position_index]
            
            # modify the list and pop name in the current cell position object
            cell_type_position_obj._all_positions[0] = layer_cell_type_position_array
            cell_type_position_obj._all_pop_names[0] = cell_type_pop_name
            
            # remove if any the other position and pop names
            cell_type_position_obj._all_positions = cell_type_position_obj._all_positions[:1] 
            cell_type_position_obj._all_pop_names = cell_type_position_obj._all_pop_names[:1] 
            
            # make sure to iterate to not pick the same coordinates twice
            last_index += desired_cell_Type_count
            
            #Select N model files according to the desired method specified in layer prop dict
            model_selection_method = population_in_layer_dict['Model_selection']
            if model_selection_method == "Uniform":
                Feature_of_importance = None
            else:
                Feature_of_importance = population_in_layer_dict['Feature']
                
            layer_abrev = layer.replace("L","")
            layer_abrev = str(layer.split("_")[1])

            model_selection_list, model_selection_df= select_models(model_features_table, 
                                                                    desired_cell_Type_count, 
                                                                    model_selection_method, 
                                                                    feature_distribution_table, 
                                                                    Feature_of_importance, 
                                                                    population, 
                                                                    layer_abrev,
                                                                    do_plot=False)

               
            layer_cell_type_node_dict = {'pop_name' : cell_type_pop_name, # "L_4_PV"
                                         'N' : desired_cell_Type_count, # 300
                                         'ei':ei, # "i"
                                         'layer':layer,  # "L_4"
                                         'cell': population, # "PV"
                                         'position_obj' : cell_type_position_obj,
                                         'model_list' : model_selection_list}
            
            layer_node_dict[population] = layer_cell_type_node_dict
        node_dict['N_tot'] = N_tot
        node_dict[layer] =  layer_node_dict  
        
    return node_dict


def function_build_nodes(node_dict, network_dir, network_name):
    
    
    net = NetworkBuilder(network_name)

    
    for current_Layer, current_layer_dict in node_dict.items():
        if current_Layer == 'N_tot':
            continue
        
        for current_population, current_pop_dict in current_layer_dict.items():
            if current_population=='Layer' or current_population == 'N':
                continue
            position_obj = current_pop_dict['position_obj']
            
            net.add_nodes(N=current_pop_dict['N'],  # Create a population of 80 neurons
                          positions = position_obj._all_positions[0],
                          
                          pop_name=current_pop_dict['pop_name'],
                          layer = current_pop_dict['pop_name'].split('_')[1],
        
                          ei= current_pop_dict['ei'],  # optional parameters
                          model_type='point_process',  # Tells the simulator to use point-based neurons
                          model_template='nest:glif_lif_asc_psc',  # tells the simulator to use NEST iaf_psc_alpha models
                          dynamics_params=current_pop_dict['model_list'] # File containing iaf_psc_alpha mdoel parameters
                         )
    
    net.build()
    net.save_nodes(output_dir=network_dir)
    

    
#%%


def add_nodes_V1_in_nrrd (dict_path,factor) : 
	"""
	Place neurons in a nrrd volume 
	--------
	dict_path : str
		Path to the dictionnary containing all necessary information to define and place the neurons within a nrrd volume (check ../Additional_data/dict_v1_nodes.json to see the dictionnary structure)
	factor : float : 
		Defines the % of the density kept to build the nodes. Has to be inferior to 1
	Returns
	-------
	nets : list
		list of bmtk object from which we will be able to build and save the nodes
	dataframes : list 
		List of dataframes containing for each the coordinates and cell type names
	"""

	v1_info=json.load(open(dict_path,'r'))
	nets=[]
	dataframes=[]
	for layer in v1_info.keys() : 
		location=v1_info[layer]["layer_name"]
		layer_factor_name=layer+"_factor_"+str(factor)
		layer_location=CellLocations(layer_factor_name)
		layer_location.dmin=v1_info[layer]["dmin"]
		layer_location.CCF_orientation=True
		maximum_density=v1_info[layer]["density"]*factor
		path_nrrd="../Components/nrrd_files/"+v1_info[layer]["nrrd_file"]
		pop_names=[]
		proportions=[]
		for subtype in v1_info[layer]["nodes_descriptions"] : 
			pop_names.append(subtype["pop_name"])
			proportions.append(subtype["proportion"])
		layer_location.add_positions_nrrd(path_nrrd,maximum_density,pop_names=pop_names,partitions=proportions,method='prog',verbose=True)
		
		
		net=NetworkBuilder(layer_factor_name)
		for i,subtype in enumerate(v1_info[layer]["nodes_descriptions"]) : 
			getattr(layer_location,pop_names[i]).tuning_angles = np.linspace(0.0,360.0,getattr(layer_location,pop_names[i]).N)	
			net.add_nodes (N=getattr(layer_location,pop_names[i]).N,
				       positions = getattr(layer_location,pop_names[i]).positions,
				       pop_name=pop_names[i],
				       model_type=subtype["model_type"],
				       model_template=subtype["model_template"],
				       dynamics_params=subtype["dynamics_params"],
				       ei=subtype["ei"],
				       location=location,
				       tuning_angle=getattr(layer_location,pop_names[i]).tuning_angles
			)
		nets.append(net)
		df_layer=get_population_loc_df(layer_location)
		dataframes.append(df_layer)

	return(nets,dataframes)

def node_interval (cell_type_name,df_csv,list_id_type) :
    """
    For a given neuron type, it returns the neuron id interval within the list of all neuron ids present in the node csv file (SONATA)
    -----
    :param cell_type_name: str
        cell type name (has to match one of the pop_name in the node csv file (SONATA)
    :param df_csv : dataframe
        the node csv file (SONATA) in a dataframe format
    :param list_id_type: list
        the node_type_id list from the node h5 file (SONATA)
    -----
    :return: list
        list of 2 values : the first and the last neuron id of the given cell type in the the node_type_id list from the node h5 file (SONATA)
    """
    df_1=df_csv.loc[df_csv["model_name"]==cell_type_name]
    node_id=int(list(df_1["node_type_id"])[0])
    pre_start=list_id_type.index(node_id)
    list_id_type.reverse()
    pre_end= len(list_id_type)-list_id_type.index(node_id)
    list_id_type.reverse()
    return([pre_start,pre_end])

def lgn_node_from_SONATA (csv_path,h5_file) : 
	"""
	Create a node/net object from an already existing lgn node SONATA file
	"""
	net=NetworkBuilder('lgn')
	h5=h5py.File(h5_file)
	csv_file=pd.read_csv(csv_path,sep=" ")
	list_split=csv_path.split("/")
	object_name=list_split[-1].split("_")
	list_id_type=list(h5["nodes"][object_name[0]]["node_type_id"][()])
	for name,subtype,model_type,model_template,dynamic_params,non_dom_params,sf_sep in zip(list(csv_file["model_name"]),list(csv_file["subtype"]),list(csv_file["model_type"]),list(csv_file["model_template"]),list(csv_file["dynamics_params"]),list(csv_file["non_dom_params"]),list(csv_file["sf_sep"])) :
        	interval=node_interval(name,csv_file,list_id_type)
        	x = list(h5["nodes"][object_name[0]]["0"]["x"][()])[interval[0]:interval[1]]
        	y = list(h5["nodes"][object_name[0]]["0"]["y"][()])[interval[0]:interval[1]]
        	tuning_angle=list(h5["nodes"][object_name[0]]["0"]["tuning_angle"][()])[interval[0]:interval[1]]
        	spatial_size=list(h5["nodes"][object_name[0]]["0"]["spatial_size"][()])[interval[0]:interval[1]]
        	net.add_nodes(
        		N=len(tuning_angle),  # n_cells,
        		model_name=name,
        		location='LGN',
        		subtype=subtype,
        		model_type=model_type,
        		model_template=model_template,
        		dynamics_params=dynamic_params,
        		non_dom_params=None if pd.isnull(non_dom_params) else non_dom_params,
        		x=x,
        		y=y,
        		tuning_angle=tuning_angle,
        		spatial_size=spatial_size,
        		sf_sep=None if pd.isnull(sf_sep) else sf_sep
    		)
	return(net)


#FUNCTION TO PLACE THE LGN CELLS ACCORDING TO A GRID
	
def is_prime(n):
  for i in range(2,n):
    if (n%i) == 0  :
      return False
  return True

def save_set_placement (liste_x,liste_y,x_grids,y_grids,liste_names,name_output_file) :
    df = pd.DataFrame()
    for ecart, names in zip(np.arange(0, len(liste_x), x_grids * y_grids), liste_names):
        ecart_x = liste_x[ecart:ecart + (x_grids * y_grids)]
        ecart_y = liste_y[ecart:ecart + (x_grids * y_grids)]
        names_x = names + "_x"
        names_y = names + "_y"
        df[names_x] = ecart_x
        df[names_y] = ecart_y
    df.to_csv(name_output_file+".csv")

def set_placement (x_grids,y_grids,lgn_models,field_size,lgn_fraction):
    n_blocks = x_grids * y_grids
    tile_width = field_size[0] / x_grids 
    tile_height = field_size[1] / y_grids 
    all_x = []
    all_y = []
    for model_name, params in lgn_models.items():
        n_cells = int(params['N']*lgn_fraction)
        per_block = n_cells / n_blocks
        prime = is_prime(int(per_block))
        xs = []
        ys = []
        if prime == True:
            per_block_bis = per_block + 1
            row = 2
            col = int(per_block_bis / 2)
            tile_width_bis = tile_width / row
            tile_height_bis = tile_height / col
            for r in range(row):
                for c in range(col):
                    xs.append(np.random.uniform(r * tile_width_bis, (r + 1) * tile_width_bis))
                    ys.append(np.random.uniform(c * tile_height_bis, (c + 1) * tile_height_bis))
            xs = xs[:-1]
            ys = ys[:-1]

        else:
            if per_block % 2 != 0:
                row = 3
            else:
                row = 2
            col = int(per_block / row)
            tile_width_bis = tile_width / row
            tile_height_bis = tile_height / col
            for r in range(row):
                for c in range(col):
                    xs.append(np.random.uniform(r * tile_width_bis, (r + 1) * tile_width_bis))
                    ys.append(np.random.uniform(c * tile_height_bis, (c + 1) * tile_height_bis))
        all_x.append(xs)
        all_y.append(ys)
    temp_x = []
    temp_y = []
    for box_x, box_y in zip(all_x, all_y):
        for i in np.arange(0, field_size[0], field_size[0] / x_grids): #16
            for j in np.arange(0, field_size[1], field_size[1] / y_grids): #12
                temp_x.append(list(np.array(box_x) + (i)))
                temp_y.append(list(np.array(box_y) + (j)))
    return(temp_x,temp_y)

def plot_grid_set (list_x,list_y,n_blocks,colors,x_grids,y_grids,azimuth_range,elevation_range) :
    for i, c in zip(np.arange(0, len(list_x), n_blocks), colors):
        ecart_x = list_x[i:i + n_blocks] #150
        ecart_y = list_y[i:i + n_blocks]


        for j, z in zip(ecart_x, ecart_y):

            plt.scatter(j, z, c=c, s=10)

    for i in np.arange(0, 127, elevation_range/y_grids): #0,132,12
        plt.axhline(i, color="black")
    for j in np.arange(0, 250, azimuth_range/x_grids): #0,256,16
        plt.axvline(j, color="black")
    plt.show()

def split_list(list_str) :
    """
    split a list of str to only keep the str numbers and store it in a list (complementary to convert_df function)
    -----
    :param list_str: list
        list of str (which has the "[","]" and "," symbols)
    -----
    :return: list
        a list of int
    """
    liste = []
    if len(list_str)==2 :
        return(liste)
    b = list_str.split("[")
    c = b[1].split("]")
    d = c[0].split(", ")
    liste=[float(i) for i in d]
    return (liste)

def convert_df (df) : #convert str lists to lists
    """
    If the lists within the dataframe are considered as str, it converts these str into a list of int (complementary to load_grid_csv function)
    :param df: dataframe
        dataframe containing str instead of list
    -----
    :return: dataframe
        dataframe containing lists (instead of str)
    """
    df=df.drop(list(df.columns)[0],axis=1)
    for i in np.arange(df.shape[0]):
        for j in list(df.columns):
            box = df.loc[i][j]
            b = split_list(box)
            df.loc[i, str(j)] = b
    return(df)

def generate_liste(name,df) :
    type_x = []
    type_y = []
    for i, j in zip(list(df[name+"_x"]), list(df[name+"_y"])):
        type_x.extend(i)
        type_y.extend(j)
    return(type_x,type_y)

def bmtk_nodes_set(lgn_models,df,out_dir,filename,lgn_fraction) : 
    set_placement = NetworkBuilder(filename)
    for model_name, params in lgn_models.items():
    	
        n_cells = int(params['N']*lgn_fraction)
        x, y = list(generate_liste(model_name, df)[0]), list(generate_liste(model_name, df)[1])
        spatial_size= params["size_range"]
        set_placement.add_nodes(
            N=len(x), 
            model_name=model_name,
            location='LGN',
            subtype=params['subtype'],
            model_type=params['model_type'],
            model_template=params['model_template'],
            dynamics_params=params['dynamics_params'],
            non_dom_params=params.get('non_dom_params', None),
            x=x,
            y=y,
            tuning_angle=[np.NaN] * len(x) if params['tuning_angle'] else np.linspace(0.0, 360.0, len(x)),
            spatial_size= [params["size_range"]]*len(x),
            sf_sep=params.get('sf_sep', None)
        )
    set_placement.build()
    set_placement.save_nodes(output_dir=out_dir)


def add_nodes_bkg(output_dir) : 
	bkg = NetworkBuilder("BKG")
	bkg.add_nodes(
	N=1, 
	pop_name='SG_001', 
   	ei='e', 
    	location='BKG',
    	model_type='virtual'
	)
	bkg.save_nodes(output_dir=output_dir)

"""
if __name__ == '__main__':
	add_nodes_bkg ("../Networks/nodes/bkg")
"""

"""
if __name__ == '__main__':
	dict_path="../Additional_data/lgn_dict.json"
	lgn_models=json.load(open(dict_path,"r"))
	lgn_fraction=0.1
	x_grids, y_grids = 5,3 #15 , 10 ; 5,3 for 10% because we want len(x_grids)*len(y_grids) =<  smallest subtype neuron number ; we want at least 1 neuron per box 
	field_size = (240.0, 120.0)	colors=["darkred","red","peru","orange","gold","olive","lime","mediumturquoise","cyan","dodgerblue","blue","darkorchid","magenta","pink"]
	liste_names=["sON_TF1","sON_TF2","sON_TF4","sON_TF8","sOFF_TF1","sOFF_TF2","sOFF_TF4","sOFF_TF8","sOFF_TF15","tOFF_TF4","tOFF_TF8","tOFF_TF15","sONsOFF_001","sONtOFF_001"]
	out_dir="../Networks/lgn/10percent"
	out_dir_fixed_coord=out_dir+"/lgn_fixed_coord"
	liste_x,liste_y=set_placement(x_grids,y_grids,lgn_models,field_size,lgn_fraction)
	plot_grid_set(liste_x,liste_y,x_grids*y_grids,colors,x_grids,y_grids,field_size[0],field_size[1])
	save_set_placement(liste_x,liste_y,x_grids,y_grids,liste_names,out_dir_fixed_coord)
	output=pd.read_csv(out_dir_fixed_coord+".csv")
	df=convert_df(output)
	bmtk_nodes_set(lgn_models,df,out_dir,"lgn",lgn_fraction)
"""

