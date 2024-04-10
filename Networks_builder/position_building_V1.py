from bmtk.builder.auxi.node_params import positions_columnar,CellLocations
import pandas as pd
from bmtk.builder.networks import NetworkBuilder
import plotly.express as px
import math
import numpy as np
import json



#%%

def create_column(param_dict, layer_list, do_plot = False):
    
    Layer_width_dict = {"L1" : 157.0,# in um
                        "L2/3" : 575.0, 
                        "L4" : 900.0,
                        "L5" : 1411.0, 
                        "L6" : 1973}
    
    Layer_center_dict = {"L1" : -78.5, #in um
                         "L2/3" : -444.5, 
                         "L4" : -1182.0,
                         "L5" : -2334.5,
                         "L6" : -4026.5}
    
    Layer_density_dict = {"L1" : 43_199.7, #from BBP, in /mm3
                         "L2/3" : 108_787.9, 
                         "L4" : 130_640.0,
                         "L5" : 105_268.4,
                         "L6" : 98_393.2}
    
    
    
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
    
    fig = px.scatter_3d(position_df, x='x', y='y', z='z',color='Cell_type')
    fig.update_traces(marker_size=4)
    fig.show('browser')
    
    return position_df

def get_population_loc_df(cell_location_obj):
    position_df=pd.DataFrame(columns=['x','y','z','Cell_type'])
    
    for pop_name,pop_loc in zip(cell_location_obj._all_pop_names,cell_location_obj._all_positions):
        
        pop_loc_df=pd.DataFrame(pop_loc)
        pop_loc_df.columns = ['x','y','z']
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
    
#%%

def add_nodes_V1_in_nrrd (dict_path,factor) : 
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
		df_layer=get_population_loc_df(layer_location)
		dataframes.append(df_layer)
		net=NetworkBuilder(layer_factor_name)
		for i,subtype in enumerate(v1_info[layer]["nodes_descriptions"]) : 
		
			net.add_nodes (N=getattr(layer_location,pop_names[i]).N,
				       positions = getattr(layer_location,pop_names[i]).positions,
				       pop_name=pop_names[i],
				       model_type=subtype["model_type"],
				       model_template=subtype["model_template"],
				       dynamics_params=subtype["dynamics_params"],
				       ei=subtype["ei"],
				       location=location,
				       tuning_angle=np.linspace(0.0,360.0,getattr(layer_location,pop_names[i]).N)
			)
		nets.append(net)
	return(nets,dataframes)


if __name__ == '__main__':
	dict_path="../Additional_data/dict_v1_nodes.json"
	net_layers,dataframes=add_nodes_V1_in_nrrd (dict_path,1) 
	plot_population(dataframes[0]) #plot l1 neurons
	#net_layers[0].build() #build l1 nodes
	#net_layers[0].save("../Networks/nodes") #save l1 nodes
	













