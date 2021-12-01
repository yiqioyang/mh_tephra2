#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:50:26 2020

@author: qingyuan
"""
# This code is used for implementing tephra2 with the Metropolis-Hastings (M-H) algorithm
	Qingyuan Yang 
		(Earth Observatory of Singapore; Asian School of the Environment, Nanyang Technological University); 
	E Bruce Pitman 
		(Department of Materials Design and Innovation and School of Engineering and Applied Sciences, University at Buffalo);
	Marcus Bursik 
		(Department of Geology, University at Buffalo);
	Susanna F Jenkins 
		(Earth Observatory of Singapore; Asian School of the Environment, Nanyang Technological University).

# Contact: qingyuan.yang@ntu.edu.sg
# License: CC BY-SA 4.0 
# Please feel free to contact us if you have any questions, comments or suggestions. 
# Last update: 21st, June, 2020


#Import packages
import sys
import os
import csv
import numpy
from numpy import *
import subprocess
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import uniform


# Set the working directory to where this code is stored        
os.chdir('/where/this/codes/is/stored/')							
                                
# Read field observations
observation = numpy.transpose(numpy.genfromtxt("../where/field/observations/are/stored/observation.txt", delimiter = " "))				

# Set up the likelihood scale
# The likelihood function is a lognormal distribution of (observation/tephra2_prediction), which 
# centers at 0 (log(1) = 0), and the value below specifies the standard deviation of the 
# lognormal distribution
likelihood_scale = 0.1

# Set up the number of runs 
runs = 10000

# Set up the check point. For example, if the variable below is 
# set to be 100, then every 100 runs, the variables to be estimated 
# will be plotted pairwise.
# If users have a lot of variables to estimate, we recommend you to 
# choose a smaller value for it in your experimental runs, and switch 
# to a larger value in the real run. This is because plotting many 
# biplots with many points could cost a lot of unnecessary pauses during the running. 
check_snapshot = 9999

############################## START OF THE NO-NEED-TO-CHANGE BLOCK ############################################################
# Set up variable names
var_nam = 'h,m,alpha,beta,max_gs,min_gs,med_gs,sd_gs,ventx,venty,ventz,edy_const,dif_coef,ftthreshld,lith_rou,pum_rou,col_steps,part_steps,plum_model,wdir,wmax_v,vmax_elev,wupper_elev'

# Set up elevations where wind speed and direction will be specified
wind_elevation = numpy.array(range(1,41))*1000

# Read the required inputs for the eruption source parameters to run the
# M-H algorithm.
esp_ini = numpy.genfromtxt("../input/esp_input.csv", 
                 delimiter = ",", skip_header=1, 
                 dtype = ("|S15", "float","|S15", "float", "float", "float"))
input_plume = numpy.transpose(esp_ini['f1'])
prior_type_plume = numpy.transpose(esp_ini['f2'])
prior_para = numpy.transpose([esp_ini['f3'],esp_ini['f4']])
draw_scale = numpy.transpose(esp_ini['f5'])

# Read the required inputs for the wind profile to run the M-H algorithm.
wind_ini = numpy.genfromtxt("../input/wind_input.csv", 
                 delimiter = ",", skip_header=1, 
                 dtype = ("|S15", "float","|S15", "float", "float", "float"))
input_wind = numpy.transpose(wind_ini['f1'])            
prior_type_wind = numpy.transpose(wind_ini['f2'])
prior_para_wind = numpy.transpose([wind_ini['f3'],wind_ini['f4']])
draw_scale_wind = numpy.transpose(wind_ini['f5'])

# Source the functions required to run the Metropolis-Hastings algorithm, which should be
# in the current path.                           
runfile("./mcmc_1.7_functions_noted.py")

# Run the algorithm and collect the output
output  = proposal_function(input_plume, input_wind, 
                          prior_type_plume, prior_type_wind,
                      draw_scale,draw_scale_wind,
                      prior_para,prior_para_wind,
                      wind_elevation,runs,likelihood_scale,observation,check_snapshot)

sample_chain = output[0]          # <- Sampled chain, MAIN RESULTS
poster_chain = output[1]          # <- Log(posterior density) for each drawn sample
accept = output[2]                # <- The number of accepted draws
prior_array = output[3]           # <- Log(prior) for each proposed sample
likeli_array = output[4]          # <- Log(likelihood) for each proposed sample
############################## END OF THE NO-NEED-TO-CHANGE BLOCK ############################################################

os.chdir('../where/you/want/to/store/the/output/')              # Change the directory to where you want to save the results

############################## START OF THE NO-NEED-TO-CHANGE BLOCK ############################################################
# Save the output into separate files. 
numpy.savetxt("sample_chain.csv", sample_chain, delimiter=',', header=var_nam)
numpy.savetxt("poster_chain.csv", poster_chain, delimiter=',', header='post')
numpy.savetxt("prior_array.csv", prior_array, delimiter=',', header='prior')
numpy.savetxt("number_of_acceptances.txt", [accept], delimiter=',')
numpy.savetxt("likeli_array.csv", likeli_array, delimiter=',', header='likelihood')
############################## END OF THE NO-NEED-TO-CHANGE BLOCK ############################################################
