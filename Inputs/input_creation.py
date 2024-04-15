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

"""
if __name__ == '__main__':
	hz=1000
	frequency_input (hz,[1,2],"inputs")
	h=h5py.File('inputs/spikes_1000Hz_1s.h5')
	times = list(h['spikes']['SpikeGenerator']['timestamps'][()])
	print('total nber of spikes over all duration : ', len(times))
"""	










