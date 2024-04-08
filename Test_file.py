#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:23:05 2024

@author: julienballbe
"""

import pandas as pd
import numpy as np
import math
import json
import matplotlib.pyplot as plt
import h5py
import glob
#import nest

from bmtk.analyzer.spike_trains import plot_rates_boxplot, plot_rates, plot_raster
from bmtk.utils import sonata
from bmtk.utils.reports import SpikeTrains
from bmtk.builder.auxi.node_params import positions_columnar
from bmtk.builder.auxi.node_params import CellLocations
#import plotly.express as px
import math

print("test Margaux bis")
