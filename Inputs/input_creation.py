import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime
from datetime import timedelta
from scipy.special import erf
import h5py
import time
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

def frequency_input (hz,stimulation_times,folder_path) : 
	
	psg=PoissonSpikeGenerator(
	population='SpikeGenerator',
	output_units='ms'
	)	
	psg.add(
	node_ids=0,
	firing_rate=float(hz),
	times=(stimulation_times[0],stimulation_times[1]),
	)
	duration=stimulation_times[1]-stimulation_times[0]
	psg.to_sonata(folder_path+'/spikes_'+str(hz)+'Hz_'+str(round(duration))+'s.h5')

def pA_input (amplitude,delay,duration,population,pop_name,config_file) : 
	"""
	Change the config_file to adapt to the given pA input
	"""
	with open (config_file,'r') as f:
		data=json.load(f)
		data["inputs"]["current_clamp"]["node_set"]["population"]=population
		data["inputs"]["current_clamp"]["node_set"]["pop_name"]=pop_name
		data["inputs"]["current_clamp"]["amp"]=float(amplitude)
		data["inputs"]["current_clamp"]["delay"]=float(delay)
		data["inputs"]["current_clamp"]["delay"]=duration
	os.remove(filename)
	with open(filename,'w') as f : 
		json.dump(data,f,indent=4)


"""
if __name__ == '__main__':
	hz=1000
	frequency_input (hz,[1,2],"inputs")
	h=h5py.File('inputs/spikes_1000Hz_1s.h5')
	times = list(h['spikes']['SpikeGenerator']['timestamps'][()])
	print('total nber of spikes over all duration : ', len(times))
"""	










