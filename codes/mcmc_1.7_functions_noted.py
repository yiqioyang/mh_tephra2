#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Created on Thu Jun 28 15:05:09 2018


# This code includes all required functions to run the dMetropolis-Hastings algorithm 
# coupled with the volcanic ash transport model tephra2. 
# Users do not need to modify anything here. 
# *** We encourage users to read notations within this script, because they provide useful
# information about this algorithm. 
# Authors:
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
# Last update: 21th, June, 2020
#=========================================================================================================
# Import necessary packages
import sys
import os
import csv
import numpy
import time
from numpy import *
import subprocess
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import uniform
import matplotlib.pyplot as plt
from scipy.stats import beta


def changing_variable(input_plume, wind_speed,wind_n0, elevation):
    # This function is used to update the initial conditions (file name: 'config.txt') and
    # wind conditions (filename: 'wind.txt') for running tephra2.
    # Input:
    #   input_plume:            a vector file of 17 components
    #       input_plume[0]:         Column height in km;
    #       input_plume[1]:         Log(total eruption mass) in kg;                      # It is the natural log
    #       input_plume[2]:	    Alpha value of the plume mass distribution;		        # New feature; See tephra2 websites for more detailed description
    #       input_plume[3]:         Beta value for the plume mass distribution;          # New feature; See tephra2 websites for more detailed description
    #       input_plume[4,5,6,7]:   Max and min of grain size distribution
    #                               and median and std of grain size distribution;
    #       input_plume[8,9,10]:     Vent coordinates and elevation;
    #       input_plume[11,12]:      Diffusion coefficient & eddy constant;               # See https://vhub.org/resources/756/download/Tephra2_Users_Manual.pdf for more details
    #       input_plume[13]:        Fall time threshold;                                 # ...
    #       input_plume[14,15]:     Lithic and pumice density;                           # ...
    #       input_plume[16]:        Column steps;                                        # ...
    #       input_plume[17]:        The number of particle bins for the GSD;             # New feature. See tephra2 websites for more detailed description
    #       input_plume[18]:        Plume model                                          # New feature. See tephra2 websites for more detailed description
    #
    #   wind_speed:             output from function "process_wind_temp"                 # *** Currently the "wind_speed" is assumed to be relatively simple; 
    #                           it specifies the wind speed at each elevation (variable name: "elevation").
    
    #   wind_n0:                wind direction.                                         # Currently assumed to be constant with elevation,
    #                                                                                   # but its value can be modified in each iteration. 
    #   elevation:              elevations that the wind speed will be specified. 
    #
    # Output:
    #   This function has no output.
    #   It only modifies files "config.txt" and "wind.txt", which are the input required to run tephra2.
    ###
    # Start to modify input file "config.txt" (initial conditions of tephra2).
    fid = open('../running_tephra2/config.txt', 'r')   
    s = fid.readlines()
    s[0] = 'PLUME_HEIGHT %f' % float(input_plume[0]) + '\n'
    s[1] = 'ERUPTION_MASS %f' % float(math.exp(input_plume[1])) + '\n'
    s[2] = 'ALPHA %f' %float(input_plume[2]) + '\n'					# New version feature
    s[3] = 'BETA %f' %float(input_plume[3]) + '\n'					# New version feature
    s[4] = 'MAX_GRAINSIZE %f' % float(input_plume[4]) + '\n'
    s[5] = 'MIN_GRAINSIZE %f' % float(input_plume[5]) + '\n'
    s[6] = 'MEDIAN_GRAINSIZE %f' % float(input_plume[6]) + '\n'
    s[7] = 'STD_GRAINSIZE %f' % float(input_plume[7]) + '\n'
    s[8] = 'VENT_EASTING %f' % float(input_plume[8]) +'\n'
    s[9] = 'VENT_NORTHING %f '% float(input_plume[9])+'\n'
    s[10] = 'VENT_ELEVATION %f' % float(input_plume[10])+'\n'
    s[11] = 'EDDY_CONST %f' % float(input_plume[11])+'\n'
    s[12] = 'DIFFUSION_COEFFICIENT %f' % float(input_plume[12])+'\n'
    s[13] = 'FALL_TIME_THRESHOLD %f' % float(input_plume[13])+'\n'
    s[14] = 'LITHIC_DENSITY %f' % float(input_plume[14])+'\n'
    s[15] = 'PUMICE_DENSITY %f' % float(input_plume[15])+'\n'
    s[16] = 'COL_STEPS %f' % float(input_plume[16])+'\n'
    s[17] = 'PART_STEPS %f' % float(input_plume[17])+'\n'				# New version feature: the number of bins for a single grain
    s[18] = 'PLUME_MODEL %d' % int(input_plume[18])+'\n'				# New version feature: can be value 2
    fid.close()

    out = open('../running_tephra2/config.txt', 'w')
    for i in range(len(s)):
        out.write(s[i])
    out.close()
    # Finished modifying file "config.txt". 
##########################################
    # Start to modify input file "wind.txt" (wind conditions).
    wind_dir = wind_n0[0]
    fid = open('../running_tephra2/wind.txt', 'r')
    s = fid.readlines()
    for i in range(len(elevation)):
        s[i]= '%f %f %f' %(float(elevation[i]),float(wind_speed[i]), float(wind_dir)) +'\n'
    
    fid.close()

    out = open('../running_tephra2/wind.txt', 'w')
    for k in range(len(s)):
        out.write(s[k])
    out.close()
    # # Finished modifying file "wind.txt".
### End of function "changing_variable"

#------------------------------------------------------------------------------------------------------------------------------
def process_wind_temp(input_wind, elevation):
    # This function is used to generate wind speed profile to run tephra2.
    # *** It is separated from function "changing_variable", because there are definitely DIFFERENT WAYS of specifying
    #     the wind speed and direction at different elevations. THE FUNCTION "process_wind_temp" IS JUST ONE SIMPLE FORM THAT IS ADOPTED TEMPORARILY.
    
    # *** Users could change the content of this function to specify more complicated wind speed profiles. 
    
    # *** Currently we use four parameters to define the wind speed profile.
    # *** We assume that the wind direction is constant with elevation, but its value can be varied in each iteration. 
    # *** We assume that the wind speed increases linearly from the ground (zero elevation, zero velocity) to a certain elevation.
    # *** and then start to decrease linearly until zero.
    # *** Based on the assumptions for tephra2, the only simplification we have made in fact is that the wind direction is constant with elevation,
    # *** but its value can be varied with iteration. 
    #
    # Input:
    #   input_wind: a vector file of four components
    #       input_wind[0]:      wind direction;
    #       input_wind[1]:      max wind speed;
    #       input_wind[2]:      the elevation with the max wind speed;
    #       input_wind[3]:      the elevation that the wind speed goes to zero;
    #   elevation: a vector file that specifies the elevations that wind speed and direction will be specified.
    # Output:
    #   wind_speed: wind speed at each elevation.
    max_speed = float(input_wind[1])
    max_speed_ele = float(input_wind[2])
    zero_windspeed_ele = float(input_wind[3])
    slope_increase = max_speed/max_speed_ele
    slope_decrease = max_speed/(zero_windspeed_ele - max_speed_ele)*-1
    
    wind_speed = zeros(size(elevation))
    
    for i in range(0, size(elevation)):
        if(elevation[i]<max_speed_ele):
            wind_speed[i] = slope_increase * elevation[i]
        else:
            wind_speed[i] = max_speed + slope_decrease * (elevation[i]-max_speed_ele)
    
    wind_speed[wind_speed<0] = repeat(0,size(wind_speed[wind_speed<0]))
    
    return (wind_speed)
### End of function "process_wind_temp"

#------------------------------------------------------------------------------------------------------------------------------
def draw_input_Gaus(input_plume, prior_type_plume, draw_scale, prior_para):
    # This function is used to draw new samples based on the proposal distribution for the initial conditions of the plume.
    # It first identifies if a certain variable for the plume (mass, column height...,) needs to be modified or not by
    # variable "prior_type_plume", and then draws the proposed samples from a Gaussian distribution centered at input_plume[i] with standard
    # deviation draw_scale[i]. Each of the initial conditions for the plume is drawn separately/independently with a for loop.
    # Input:
    #     input_plume:            Initial conditions for the plume at the jth iteration;
    #                             it should look like/has the same dimensions as the "input_plume" in function "changing_variable".
    #     prior_type_plume:       The prior type for the initial conditions for the plume;             
    #                             if certain initial conditions are assumed to be known, then it is "Fixed", meaning that we keep
    #                             it fixed, and do not draw new values for it.
    #                             See "run_mh.py" file for an example of this variable.
    #     draw_scale:             The standard deviation of the proposal function (a Gaussian distribution);
    #                             See "run_mh.py" file for an example of this variable.
    #     prior_para:             The parameters that are used to define the prior;
    #                             For example, if prior_type_plume[0] = "Gaussian", and prior_para[[0]] = [0,1]
    #                             this specifies the prior for input_plume[0] will be a Gaussian distribution centered at 0 with std being 1.
    #                        
    # Output:
    #    plume_draw:             Proposed samples drawn based on the Gaussian distribution centered at input_plume[i] and draw_scale[i];
    #                            It is a vector that has the same dimesions as variable "input_plume" with the variables that are "Fixed" unchanged.
    # *** Note:
    #    Currently the code only supports Gaussian and Uniform distributions for the prior.
    plume_draw = zeros(len(input_plume))
    for i in range(len(input_plume)):
        para_local = numpy.array(prior_para[i])
        if prior_type_plume[i] == "Gaussian":
            plume_draw[i] = numpy.random.normal(loc = input_plume[i], scale = draw_scale[i], size = 1)	
            while norm.pdf(plume_draw[i], loc = para_local[0], scale = para_local[1]) == 0:
                plume_draw[i] = numpy.random.normal(loc = input_plume[i], scale = draw_scale[i], size = 1)     
                print "try resample plume"
                
        if prior_type_plume[i] == "Uniform":
            plume_draw[i] = numpy.random.normal(loc = input_plume[i], scale = draw_scale[i], size = 1)	
            while plume_draw[i] < para_local[0] or plume_draw[i] > para_local[1]:
                plume_draw[i] = numpy.random.normal(loc = input_plume[i], scale = draw_scale[i], size = 1)
                print "try resample plume"
                
        if prior_type_plume[i] == "Fixed":
            plume_draw[i] = input_plume[i]
    return plume_draw
#   End of function "draw_input_Gaus"

#------------------------------------------------------------------------------------------------------------------------------
def draw_wind_Gaus(input_wind, prior_type_wind, draw_scale_wind, prior_para_wind):
    # This function is used to draw new samples based on the proposal distribution for the wind conditions.
    # It first identifies if a certain variable for the wind conditions needs to be modified or not by
    # variable "prior_type_wind", and then draws the proposed samples from a Gaussian distribution centered at input_wind[i] with standard
    # deviation draw_scale_wind[i]. Each parameter for the wind conditions is drawn separately/independently with a for loop.
    # Input:
    #     input_wind:             Wind conditions at the jth iteration;                                              
    #                             it should look like/has the same dimensions as the "input_wind" in function "process_wind_temp".
    #     prior_type_wind:        The prior type for the wind conditions;             
    #                             if certain initial conditions are assumed to be known, then it is "Fixed", meaning that we keep
    #                             it fixed, and do not draw new values for it.
    #                             See "run_mh.py" file for an example of this variable.
    #     draw_scale_wind:        The standard deviation of the proposal function (a Gaussian distribution);
    #                             See "run_mh.py" file for an example of this variable.
    #     prior_para_wind:        The parameters that are used to define the prior.
    #                             For example, if prior_type_wind[0] = "Gaussian", and prior_para_wind[[0]] = [0,1],
    #                             that specify the prior for the input_wind[0] to be a Gaussian distribution centered at 0 with std being 1;
    #                             
    # Output:
    #     wind_draw:              Proposed samples drawn based on the Gaussian distribution centered at input_wind[i] and draw_scale_wind[i];
    #                             It is a vector that has the same dimesions as variable "input_plume"
    # *** Note:
    #     Currently the code only supports Gaussian and Uniform distributions for the prior.

    wind_draw = zeros(len(input_wind))
    for i in range(len(input_wind)):
        para_local = numpy.array(prior_para_wind[i])
        if prior_type_wind[i] == "Gaussian":
            wind_draw[i] = numpy.random.normal(loc = input_wind[i], scale = draw_scale_wind[i], size = 1)
            if i == 1:
                while wind_draw[i] <0:
                    wind_draw[i] = numpy.random.normal(loc = input_wind[i], scale = draw_scale_wind[i], size = 1)
                    print "drawn wind speed below zero"
                    
        if prior_type_wind[i] == "Uniform":
            wind_draw[i] = numpy.random.normal(loc = input_wind[i], scale = draw_scale_wind[i], size = 1)
            while wind_draw[i] < para_local[0] or wind_draw[i] > para_local[1]:
                wind_draw[i] = numpy.random.normal(loc = input_wind[i], scale = draw_scale_wind[i], size = 1)
                print "try resample wind"
                
            if i == 1:
                while wind_draw[i] <0:
                    wind_draw[i] = numpy.random.normal(loc = input_wind[i], scale = draw_scale_wind[i], size = 1)
                    print "drawn wind speed below zero"
                
        if prior_type_wind[i] == "Fixed":
            wind_draw[i] = input_wind[i]
    return wind_draw
#   End of function "draw_wind_Gaus"

#------------------------------------------------------------------------------------------------------------------------------
def likelihood_function(prediction, observation, likelihood_ratio):
    # This function calculates value of the likelihood function for each pair of observation and simulation;
    # Input:
    #       prediction:             Simulated values from tephra2.
    #       observation:            Observations from the field;
    #                               In our example, the ``observations'' are generated from tephra2.
    #       likelihood_ratio:       A constant that characterizes the scale of the likelihood function.
    # Output:
    #       likelihood_array:       A vector that contains the LOG-TRANSFORMED value of the likelihood function for each pair
    #                               of observation and simulation.
    # *** It should be noted that different forms of likelihood function can be specified.
    # *** The one we are using scales with the magnitude of observationm and is a 
    # *** Gaussian for the value log(observation/prediction) with a mean of 0, and scale of "likelihood_ratio".
    # *** We assume that each observation is made independently;
    # *** If prediction from tephra2 produces value below 0.001, it is turned to 0.001 to avoid "any value/0" problem.
    n = len(observation)                                          # Number of observations
    likelihood_array = zeros(n)
    prediction[prediction<0.001] = repeat(0.001, size(prediction[prediction<0.001]))
        # ****** We turn simulations that are thinner than 0.001 (including 0) to 0.001 to avoid negative infinity
        #        for the log of the likelihood function because log(0) goes to infinity. 
    for i in range(n):
        likelihood_array[i] = math.log(norm.pdf(math.log(observation[i]/prediction[i],10),loc = 0, scale = likelihood_ratio),10)
        # ****** Please note the FORMAT of the likelihood function. It SCALES with the observation. 
    return likelihood_array
#   End of function "likelihood_function"
#------------------------------------------------------------------------------------------------------------------------------

def prior_function_plume(prior_type_plume, input_plume, prior_para):		
    # This function calculates the prior for the initial conditions of the plume
    # Input:
    #       prior_type_plume:   The prior types for the initial conditions of the plume;
    #                           Currently the code only supports Uniform and Gaussian distributions;
    #                           If some initial conditions on the plume are known to be fixed (e.g., we know exactly the vent location),
    #                           we just keep it "Fixed", and not calculate the prior;
    #                           See "run_mh.py" for an example of this variable.
    #       input_plume:        The initial conditions of the plume at the jth iteration.
    #       prior_para:         The parameters for the prior, if prior_type_plume[0] = "Gaussian", prior_para[0] = [0,1] specifies
    #                           that it follows a Gaussian distribution centered at 0 with a std of 1.
    #                           If prior_type_plume[0] = "Uniform", prior_para[0] = [0,1], the [0,1] specifies the lower and upper
    #                           bounds for the uniform distribution.
    #                           See "run_mh.py" for an example of this variable.
    # Output:
    #       priors:             The log-transformed prior (base:10) for each initial condition for the plume that needs to be 
    #                           estimated.
    n = len(input_plume)
    priors = zeros(n)
    for i in range(n):
        if prior_type_plume[i] == 'Uniform':
            para_local = numpy.array(prior_para[i])
            priors[i] = math.log(uniform.pdf(input_plume[i], loc = para_local[0], scale = para_local[1]-para_local[0]),10)
                
        if prior_type_plume[i] == 'Gaussian':
            para_local = prior_para[i]
            priors[i] = math.log(norm.pdf(input_plume[i], loc = para_local[0], scale = para_local[1]),10)

        if prior_type_plume[i] == 'Fixed':
            priors[i] = 0
            
    return (priors)
#   End of  function "prior_function_plume"
#------------------------------------------------------------------------------------------------------------------------------
def prior_function_wind(prior_type_wind, input_wind, prior_para_wind):
    # This function calculates the prior for the wind conditions.
    # Input:
    #       prior_type_wind:    The prior types for the wind conditions;
    #                           Currently the code only supports Uniform and Gaussian distributions;
    #                           If some wind conditions are known to be fixed (e.g., we know exactly the wind direction),
    #                           we just keep it "Fixed", and not calculate the prior;
    #                           See "run_mh.py" for an example of this variable.
    #       input_wind:         The wind conditions at the jth iteration.
    #       prior_para_wind:    The parameters for the prior, if prior_type_plume[0] = "Gaussian", prior_para[0] = [0,1] specifies
    #                           that it follows a Gaussian distribution centered at 0 with a std of 1.
    #                           If prior_type_plume[0] = "Uniform", prior_para[0] = [0,1], the [0,1] specifies the lower and upper
    #                           bounds for the uniform distribution.
    #                           See "run_mh.py" for an example of this variable.
    # Output:
    #       priors_wind:        The log-transformed prior (base:10) for each wind condition that needs to be 
    #                           estimated.
    
    n = len(input_wind)
    priors_wind = zeros(n)
    for i in range(n):
        if prior_type_wind[i] == 'Uniform':
            para_local = prior_para_wind[i]
            priors_wind[i] = math.log(uniform.pdf(input_wind[i], loc = para_local[0], scale = para_local[1]-para_local[0]),10)

        if prior_type_wind[i] == 'Gaussian':
            para_local = prior_para_wind[i]
            priors_wind[i] = math.log(norm.pdf(input_wind[i], loc = para_local[0], scale = para_local[1]),10)
       
        if prior_type_wind[i] == 'Fixed':
            priors_wind[i] = 0
        
    return (priors_wind)
#------------------------------------------------------------------------------------------------------------------------------
def proposal_function(input_plume, input_wind, 
                      prior_type_plume, prior_type_wind,
                      draw_scale,draw_scale_wind,
                      prior_para,prior_para_wind,
                      elevation,runs,likelihood_scale,observation,check_snapshot):          
    # This function makes use of all above functions, and implments the Metropolis-Hastings algorithm;
    # Input:
    #       All input variables except for "runs" have been introduced in previous functions;
    #       See "run_mh.py" for examples of these variables.
    #       runs:       The number of sample draws for the Metropolis-Hastings algorithm. 
    # Output:
    #       chain:              The chain that contains all samples drawn for the Metropolis-Hastings algorithm;
    #                           The initial conditions for the plume and wind conditions that are specified to be "Fixed"
    #                           are kept in the chain (elments in the corresponding column have the same value).
    #       post_chain:         Value of likelihood function*prior for each accepted/kept sample.
    #       acceptance_count:   The number of samples that are accepted during the implmentation of the Metropolis-Hastings algorithm.
    #                           Acceptance rate = acceptance_count/runs.
    #                           This can be used to evaluate the performance of the method.
    #       prior_array:        log(prior) for each draw.
    #       likeli_array:       log(likelihood function) for each draw.
    #
    #*** KEY ASSUMPTIONS for the present version of the code:
    #   *** Each initial condition has indepdent prior distribution.
    #   *** Each observation is made independently, such that the likelihood function could be summed up under log-scale.
    #   *** Our likelihood function is just one possible form of likelihood functions.
    #   *** Wind direction is constant with elevation, but its value can be changed in each iteration (we could estimate the wind direction).
    # Additional notes to this function are added below.
    #------------------------------------------
    # Beging preparation for the plot
    var_name = ['h', 'm', 'alpha', 'beta','max_gs', 'min_gs', 'med_gs', 'sd_gs', 'ventx', 'venty', 'ventz', 'edy_const', 'dif_coef', 'ftthreshld', 'lith_rou', 'pum_rou','col_steps','part_steps', 'plum_model', 'wdir', 'wmax_v', 'vmax_elev','wupper_elev']
    prior_pl = numpy.asarray(prior_type_plume)
    prior_wd = numpy.asarray(prior_type_wind)
    
    n_var = len(prior_pl[prior_pl!="Fixed"]) + len(prior_wd[prior_wd!="Fixed"])
    
    var_index_pl = numpy.where(prior_pl!="Fixed")
    var_index_wd = numpy.add(numpy.where(prior_wd!="Fixed"),19)
    var_index = numpy.append(var_index_pl,var_index_wd)
    
    #chain = zeros([runs+1,len(input_plume+input_wind)])                     # Specify an empty matrix to store the drawn samples (the sampled markov chain) from the M-H algorithm.
    chain = zeros([runs+1,len(numpy.append(input_plume,input_wind))])  
    
    post_chain = zeros(runs+1)                                              # Specify an empty vector to store the value of the likelihood function*prior for each accepted/drawn sample.
    acceptance_count = 0                                                    # Specify a variable to count the number of acceptance for the method.
    prior_array = zeros(runs+1)                                             # Specify an empty vector to store the log(prior) for each draw.
    likeli_array = zeros(runs+1)                                            # ... to store the log(likelihood function) ...  
    
    input_plume_n0 = input_plume                                            # Specify the initial value for the initial conditions for the plume.
    input_wind_n0 = input_wind                                              # Specify the initial value for wind conditions.
    #-------Begin of section A-----------------#
    wind_speed_n0 = process_wind_temp(input_wind, elevation)                # Turn our current way of specifying the wind conditions to wind speed and direction and each elevation.       
    #chain[0,] = input_plume_n0 + input_wind_n0                              # Store the initial value to the first row of the chain.
    chain[0,] = numpy.append(input_plume_n0, input_wind_n0)                              # Store the initial value to the first row of the chain.
    
    
    changing_variable(input_plume_n0, wind_speed_n0, input_wind, elevation) # Update the specified initial conditions and wind condtions required to run tephra2.
    subprocess.call(["../running_tephra2/jobs_wcf.sh"])
                                                                            # Run tephra2
    prediction = numpy.transpose(numpy.transpose(numpy.genfromtxt("../running_tephra2/a_t2_output.txt", delimiter = ' '))[[3]])
                                                                            # Collect simulations from tephra2
    likelihood_temp = likelihood_function(prediction, observation, likelihood_scale)
                                                                            # Calculate values of the likelihood function 
    prior_plume = prior_function_plume(prior_type_plume, input_plume_n0, prior_para)
                                                                            # Calculate the prior for the initial conditions of the plume
    prior_wind = prior_function_wind(prior_type_wind, input_wind_n0, prior_para_wind)
                                                                            # Calculate the prior for the wind conditions
    posterior_temp = sum(likelihood_temp) + sum(prior_plume) + sum(prior_wind)
                                                                            # Calculate the log(prior*likelihood function)(it is proportional to the posterior distribution)
    #-------End of section A-----------------#
    post_chain[0] = posterior_temp 
                                                                            # Store the log(prior*likelihood function) for the initially specified value
    input_plume_temp = input_plume_n0                                       # Set the specified initial conditions for the plume to be a temporal variable, such that we could just
                                                                            # update this temporal variable in the following iterations.
    input_wind_temp = input_wind_n0                                         # Set the specified wind conditions to be a temporal variable, such that we could just update 
                                                                            # this temporal variable in the following iterations.
                                                                            # In this way, the following iterations will always update variables *_temp to avoid confusion.
    n_variable = len(prior_type_plume[prior_type_plume != 'Fixed']) + len(prior_type_wind[prior_type_wind != 'Fixed'])
    
    white_space = "          "*n_variable                                   # Line used to monitor the implementation of the method
                                                                            # ...
    print("Propose%s | Result \n" %white_space)                             # ...
    
    prior_array[0] = sum(prior_plume) + sum(prior_wind)                     # Store the log(prior) in prior_array.
    likeli_array[0] = sum(likelihood_temp)                                  # Store the log(likelihood) in likeli_array.
    for i in range(1,(runs+1)):                                             # Start of the iteration to draw samples.

            
        input_plume_temp = draw_input_Gaus(chain[i-1,range(19)], prior_type_plume, draw_scale,prior_para)
                                                                            # Propose a new draw for the initial conditions of the plume based on the proposal distribution and
                                                                            # values from the previous iteration.
        input_wind_temp = draw_wind_Gaus(chain[i-1,range(19,23)], prior_type_wind,draw_scale_wind,prior_para_wind)
                                                                            # Propose a new draw for the wind conditions based on the proposal distribution and values from the previous
                                                                            # iteration.
    #----------Begin of section B--------------#
    # This section is identical to Section A of the code just a few lines above                                                                        
        wind_speed_temp = process_wind_temp(input_wind_temp, elevation)     #
        changing_variable(input_plume_temp, wind_speed_temp, input_wind_temp , elevation)
        subprocess.call(["../running_tephra2/jobs_wcf.sh"])				
        pred_temp = numpy.transpose(numpy.transpose(numpy.genfromtxt("../running_tephra2/a_t2_output.txt", delimiter = ' '))[[3]])	
        likelihood_temp = likelihood_function(pred_temp, observation, likelihood_scale)
        prior_plume = prior_function_plume(prior_type_plume, input_plume_temp, prior_para)   
        prior_wind = prior_function_wind(prior_type_wind, input_wind_temp, prior_para_wind)
        posterior_temp = sum(likelihood_temp) + sum(prior_plume) + sum(prior_wind)
    #----------End of section B--------------#
        post_ratio = posterior_temp - post_chain[i-1]                       # Calculate the logged-ratio of the posterior distribution of the current proposed sample (at the ith iteration) to the
                                                                            # to the posterior distribution of the (i-1)th sample.
        acceptance_ratio= math.log(numpy.random.uniform(low=0, high=1.0, size=1), 10)
                                                                            # Propose a random value with unifrom distribution from 0-1, and take its logarithm. 
        
        if post_ratio > acceptance_ratio:                                   # Compare the post_ratio to the acceptance_ratio.
                                                                            # We keep/accept the proposed sample, and set the ith row of the chain as the proposed value.
                                                                            # posterier > the posterier from the last time, we accept it
            chain[i,] = numpy.append(input_plume_temp, input_wind_temp)
            post_chain[i] = posterior_temp
            acceptance_count = acceptance_count+1
            accept_or_not = "accept"

        else:                                                               # We reject the proposed sample, and set the ith row of the chain as
                                                                            # the (i-1)th row of the chain.                                                                       
            chain[i,] = chain[i-1,]
            post_chain[i] = post_chain[i-1]
            accept_or_not = " rejct"

        prior_array[i] = sum(prior_plume) + sum(prior_wind)                 # Store the log(prior) in prior_array.
        likeli_array[i] = sum(likelihood_temp)                              # Store the log(likelihood) in likeli_array.
        
        print_proposed = numpy.append(input_plume_temp[numpy.where(numpy.transpose(prior_type_plume) != 'Fixed')],
                         input_wind_temp[numpy.where(numpy.transpose(prior_type_wind) != 'Fixed')])
                                                                            # Store the proposed values for printout during the ith iteration.
                                                                            # Line used to monitor the performance of the method during implementation. 
        print_old_chain = numpy.append(chain[i-1, numpy.where(numpy.transpose(prior_type_plume)!= 'Fixed')],
                                        chain[i-1, numpy.array(numpy.where(numpy.transpose(prior_type_wind) != 'Fixed'))+19])
         
        # Lines below are all used to monitor the performance of the method, and would not affect the results. 
        # Please feel free to modify them or comment them out when you are using this method. 
        # Modifiable lines below
        #
        if i % check_snapshot == 0:                                                     # Line used to monitor the performance of the method.we
            # Every check_snapshot iterations, the method gives a snapshot (see notations below) of the chain and 
            # what is going on at the check_snapshot*1, check_snapshot*2,.., th iteration.
            # if check_snapshot = 100, the method will show the snapshot every 100 iteratons.
            # Showing the snapshot more frequently could slow down the speed of the method (plotting >thousands of points).
            # It is suggested to change check_snapshot to a small value during experiments.
            # Once you are sure that the method could run properly, you could change check_snapshot to a greater value, and just focus on the final output.
            print(i)                                                        # Line used to monitor the performance of the method.
            plt.figure(figsize=(10,10))
            for m in range((n_var-1)):
                for n in range(m+1,n_var):
                    plt.subplot((n_var-1), (n_var-1), m*(n_var-1)+n)
                    plt.plot(chain[range(i+1),var_index[m]], chain[range(i+1),var_index[n]], 'o-', alpha = 0.05) 
                    plt.plot(chain[i,var_index[m]], chain[i,var_index[n]], 'ro', alpha = 0.5)
                    plt.xlabel(var_name[var_index[m]])
                    plt.ylabel(var_name[var_index[n]])
            plt.show()
        #                                                                   # These lines above specify that after each 
                                                                            # 10s of runs, the sampled variables (that are specified to be non-"Fixed") 
                                                                            # will be plotted pairwise.                                                                             
         #    time.sleep(5)
            for k in range(len(print_old_chain)):
                print("%8.2f" %print_old_chain[k]),;                            # Print out the (i-1)th row of the chain (values that are to be estimated).
            print("")
            for k in range(len(print_proposed)):                                # Print out the proposed values (values that are to be estimated) for the ith iteration.
                print("%8.2f" %print_proposed[k]),;
            print(" *** %s; (%6.2f)-(%6.2f) = %6.2f;  %4.2f"  %(accept_or_not, posterior_temp, post_chain[i-1], post_ratio, acceptance_ratio));
                                                                            # Print out: (1) if the proposed value is accepted or not;
                                                                            # (2) the log(likelihood function*prior) for the proposed sample (at ith iteration), 
                                                                            # (3) the log(likelihood function*prior) for the (i-1)th row of the chain,
                                                                            # (4) their difference (ratio in real scale) for comparison, and
 #       print("----------" * len(print_proposed))                          # (5) acceptance ratio.
        plt.close()
    return chain, post_chain, acceptance_count, prior_array, likeli_array


#------------------------------------------------------------------------------------------------------------------------------
