import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
import numpy as np
from PIL import Image as im
import pandas as pd
import math
import sys

def show_movie(movie_file, frames):
    """Helps visualize the movie"""
    movie_array = np.load(movie_file)
    fig, ax = plt.subplots(1, len(frames), figsize=(40, 5*len(frames)))

    for i, frame in enumerate(frames):
        ax[i].imshow(movie_array[frame, :, :], cmap='gray', vmin=-1.0, vmax=1.0)
        # ax[i].set_xticks([])
        ax[i].set_xticks([0])
        ax[i].set_xticklabels([frame],fontsize=20)

        ax[i].set_yticks([])

   # ax[5].set_xlabel('Frame #', horizontalalignment='right', fontsize=20)
    plt.subplots_adjust(wspace=0, hspace=0)

    plt.show()


def build_movie_template (azimuth_axis,elevation_axis,n_frames,simulation_duration) : #n_frames = 1000 : 1 image per ms = 1s
    image_path = "movies/gray_test.png"
    mat_output = np.zeros((n_frames, elevation_axis, azimuth_axis), dtype=float)
    mat_iter = 0
    img_png_path = image_path
    pic = im.open(img_png_path).convert('L')
    pic = pic.resize((azimuth_axis, elevation_axis))
    pic_data = np.asarray(pic)
    pic_data = pic_data.astype(dtype=float) * 2.0 / 255.0 - 1.0
    mat_output[mat_iter:mat_iter + simulation_duration, :, :] = pic_data
    mat_iter += simulation_duration
    return(mat_output)

def full_gray_movie (simulation_duration,mat_output,azimuth,elevation) : #the para have to be the same as in build_movie_template
    for i in np.arange(simulation_duration[0],simulation_duration[1],simulation_duration[2]) : #simulation_duration is (start,end,step)
        for j in np.arange(elevation):
            for k in np.arange(azimuth) :
                mat_output[i][j][k]=0.0
    return(mat_output)

def full_vertical_strips (strip_width,stimulus_duration,mat_output) : 
	for i in np.arange(stimulus_duration[0],stimulus_duration[1],stimulus_duration[2]) :
	    for k in np.arange(0,azimuth,strip_width*2): # k = start of each black strip
	        for u in np.arange(k,k+strip_width,1) : #u = each azimuth value within black strip
	            for j in np.arange(elevation) :
	                if u<azimuth :
	                    mat_output[i][j][u]=-1
	    for k in np.arange(strip_width,azimuth,strip_width*2):
	        for u in np.arange(k,k+strip_width,1) :
	            for j in np.arange(elevation) :
	                if u<azimuth :
	                    mat_output[i][j][u]=1
	return(mat_output)

def full_horizontal_strips (strip_width,stimulus_duration,mat_output) : 
	for i in np.arange(stimulus_duration[0],stimulus_duration[1],stimulus_duration[2]) :
	    for k in np.arange(0,elevation,strip_width*2): # k = start of each black strip
	        for u in np.arange(k,k+strip_width,1) : #u = each azimuth value within black strip
	            for j in np.arange(azimuth) :
	                if u<elevation :
	                    mat_output[i][u][j]=-1
	    for k in np.arange(strip_width,elevation,strip_width*2):
	        for u in np.arange(k,k+strip_width,1) :
	            for j in np.arange(azimuth) :
	                if u<elevation :
	                    mat_output[i][u][j]=1
	return(mat_output)

def vertical_gabor (radius,center,two_strips_width,stimulus_duration,mat_output) : #the mat_output in input should be mat_output_gray
	starts=np.arange(center[0]-radius,center[0]+radius,two_strips_width) #all azimuth values at the start of black strip within circle 
	liste=[] #all azimuth values within each first strip (black or white)
	for i in starts :
    		a=list(np.arange(i,i+(two_strips_width/2)))
    		liste.extend(a)
	for i in np.arange(stimulus_duration[0],stimulus_duration[1],stimulus_duration[2]) :
    		for j in np.arange(elevation):
        		for k in np.arange(azimuth) :
            			circle=np.sqrt(((k-center[0])**2)+((j-center[1])**2))
            			if circle<radius :
                			if k in liste :
                    				mat_output[i][j][k] = -1
                			else :
                    				mat_output[i][j][k]=1
            			else :
                			mat_output[i][j][k]=0.0
	return(mat_output)

def horizontal_gabor (radius,center,two_strips_width,stimulus_duration,mat_output) : #the mat_output in input should be mat_output_gray
	starts=np.arange(center[1]-radius,center[1]+radius,two_strips_width) #all azimuth values at the start of black strip within circle 
	liste=[] #all azimuth values within each first strip (black or white)
	for i in starts :
    		a=list(np.arange(i,i+(two_strips_width/2)))
    		liste.extend(a)
	for i in np.arange(stimulus_duration[0],stimulus_duration[1],stimulus_duration[2]) :
    		for j in np.arange(azimuth):
        		for k in np.arange(elevation) :
            			circle=np.sqrt(((j-center[0])**2)+((k-center[1])**2))
            			if circle<radius :
                			if k in liste :
                    				mat_output[i][k][j] = -1
                			else :
                    				mat_output[i][k][j]=1
            			else :
                			mat_output[i][k][j]=0.0
	return(mat_output)

def full_colored_screen (stimulus_duration,color,mat_output,elevation,azimuth) : #if not white, it's black
	if color =="white" :
		c=1.0
	else:
		c=-1.0
	for i in np.arange(stimulus_duration[0],stimulus_duration[1],stimulus_duration[2]) :
		for j in np.arange(elevation):
			for k in np.arange(azimuth) :
				mat_output[i][j][k] = c
	return(mat_output)

if __name__ == '__main__':
	azimuth=240
	elevation=120
	stimulus_duration=[400,600,1] # (start,end,step) has to be less than the duration of the simulation
	simulation_duration=[0,1000,1]
	mat_output_template=build_movie_template (azimuth,elevation,1000,1000) # imulation of 1000ms (simulation_duration) with 1image per ms (n_frames)
	
	#full gray movie (baseline for neurons)
	mat_output_gray=full_gray_movie (simulation_duration,mat_output_template,azimuth,elevation)
	
	#full grating screen
	strip_width=15
	#full_vertical_strips (strip_width,stimulus_duration,mat_output_gray) 
	#full_horizontal_strips (strip_width,stimulus_duration,mat_output_gray) 
	
	#gabor grating
	radius=32
	center=[48,60] #[60,60] left ; [180;60] right
	two_strips_width=10
	#vertical_gabor (radius,center,two_strips_width,stimulus_duration,mat_output_gray)
	#horizontal_gabor (radius,center,two_strips_width,stimulus_duration,mat_output_gray)
	
	#full colored screen (white or black) : 
	color="black"
	mat_output=full_colored_screen(stimulus_duration,color,mat_output_gray,elevation,azimuth)
	
	#save and plot movie
	movie_path="movies/full_black.npy"
	#np.save(movie_path, mat_output)
	#show_movie(movie_file=movie_path, frames=range(0, 1000, 100))


