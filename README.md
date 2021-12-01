==========================================================================================================================================
The files in this directory construct a working version of the Metropolis-Hastings algorithm 
coupled with volcanic ash transport model Tephra2.
It can be used to estimate the eruption source parameters of tephra fall deposits based on 
thickness/mass per unit area measurements in the field. 

This README file provides information about the algorithm, scripts and files
that are associated with it, and instructions on what to prepare to implement the
algorithm. 

*** NOTE: two versions of the algorithm exist. One is for a relatively simplified wind profile (files within this directory)
*** 	  and the other is for a complex wind profile.
*** 	  Their differences are: for the former, the wind direction is assumed to be constant with elevation,
*** 	  and the wind speed increases linearly with elevation until a certain elevation, and then decreases to zero
*** 	  with elevation;
***	  For the latter, wind speed and direction can be any reasonable values at each specified elevation. 
*** 	  Files within this path are served for the first way to specify the wind profile.

Authors: 
	Qingyuan Yang 
		(Earth Observatory of Singapore; Asian School of the Environment, Nanyang Technological University); 
	E Bruce Pitman 
		(Department of Materials Design and Innovation and School of Engineering and Applied Sciences, University at Buffalo);
	Marcus Bursik 
		(Department of Geology, University at Buffalo);
	Susanna F Jenkins 
		(Earth Observatory of Singapore; Asian School of the Environment, Nanyang Technological University).

Contact: qingyuan.yang@ntu.edu.sg
Last update: 21st, June, 2020
License: CC BY-SA 4.0
Please feel free to contact us if you have any questions, comments, or suggestions.

In the following text, we assume below that the path "./" refers to the path of "mcmc_Tephra2_simplifiedwind_v1.7".
------------------------------------------------------------------------------------------------------------------------------------------
Method and model:
	The Metropolis-Hastings algorithm is one of the most popular Markov Chain Monte Carlo methods. 
	It works under the Bayesian framework, and can be used to estimate the posterior distribution for variables of interest 
	(e.g., eruption source parameters). 
	There are a lot of resources on the theory and implementation of the Metropolis-Hastings algorithm. 
	See, for example, the book "Statistical and Computational Inverse Problems" by Jari Kaipio and Erkki Somersalo. 

	Tephra2 is a widely-used semi-analytical model which predicts the thickness/mass per 
	unit area of tephra fall deposits at specified locations.
	See https://vhub.org/resources/756/download/Tephra2_Users_Manual.pdf, 
	and https://github.com/geoscience-community-codes/Tephra2 for more details on Tephra2.
------------------------------------------------------------------------------------------------------------------------------------------
Input and output: 
	Input of the current version of this algorithm includes two files, and they should be stored at ./input/
	They are:
		esp_input.csv: it defines required specifications to run the M-H algorithm for the plume (Eruption Source Parameters;
		ESPs);
		wind_input.csv:it defines required specifications to run the M-H algorithm for the wind speed profile.
	Both files specify what are the variables to be estimated, and what variables are known;
	If they are known, their values need to be provided in these input files. 
	For the variables to be estimated, the prior type (e.g., Uniform distribution) and its associated parameters need to be 
	specified. Scales of the proposal function also need to be specified for each variable. 
	
	See "Instructions on what to prepare for to run the algorithm" below for more details.
	See files in the path ./input/ for input file examples.

	Output of the algorithm for this version includes:
		Sampled chain, MAIN RESULTS;
		They have 23 columns:
			The first 19 of them are ESPs for Tephra2;
			The other 4 are 
				Wind direction;	
				Maximum wind speed;
				Elevation with the maximum wind speed;
				Elevation where wind speed reaches zero.
		*** If some of the variables are known, i.e., "Fixed",
		*** their values will be fixed in the corresponding column.
		Log(posterior density) for each drawn sample;
		The number of accepted draws;	
		Log(prior) for each proposed sample;
		Log(likelihood) for each proposed sample.

------------------------------------------------------------------------------------------------------------------------------------------
Advantages:
We consider the presented algorithm having the following advantages:
	(1) Have the ability to quantify uncertainty
	The Metropolis-Hastings algorithm estimates the ESPs and wind conditions as 
	posterior probability distributions.
	
	(2) User-friendliness
	Users only need to  	(1) put the observation file and its corresponding sample sites in the appointed paths,  
				(2) modify input files in the path "./input/", and 
				(3) run the code "./codes/run_1.7_codes.py" to implement the algorithm.
	
	(3) Education 
	Tephra2 is an efficient and widely-used geophysical model, and the Metropolis-Hastings algorithm is 
	one of the most widely-used Markov Chain Monte Carlo methods.
	We hope that the presented algorithm could be used as an in-class example/tool to help students better 
	understand and have hands-on experience in using inversion methods solving geophysical problems. 	
------------------------------------------------------------------------------------------------------------------------------------------
File structure and description:
	Files in ./codes:
		Python scripts used to implement the M-H algorithm
		mcmc_1.7_functions_noted.py:
			Functions required to run the M-H algorithm.
			Users do not need to do anything with this code.
			Notations within this code are given, and we encourage interested users to read the notations within.
			
***		run_1.7_codes.py        *** The script that users run to implement the algorithm.*** 
***				        *** Please read notations within it to run it properly.
		node_ and plume2.dat are intermediate files produced from running Tephra2. 

	Files in ./input
		Input files required to run the algorithm;
		They are: esp_input.csv and wind_input.csv as introduced earlier.

	Files in ./output
		Where results will be/are stored.

	Files in ./running_tephra2
		Files required to run Tephra2.  
		In each iteration, the config.txt and wind.txt files will be modified by the algorithm
		to change the ESPs and wind profile.
		
	Files in ./src
		Source codes for Tephra2.
		
------------------------------------------------------------------------------------------------------------------------------------------
Instructions on what to prepare for to run the algorithm:
	Prerequisite:
		Linux, macOS, or Linux-like environment (as long as Tephra2 can be run by the command line).
	
	0a. Install python and spyder or other type of python IDE on the computer;

	0b. Download the file "mcmc_Tephra2_simplifiedwind_v1.7";
	
	0c. Download and install Tephra2 on your own computer (https://github.com/geoscience-community-codes/Tephra2).

	1. Navigate to the directory to the path of mcmc_Tephra2_simplifiedwind_v1.7;

	2. Delete all files in ./src, and put your installed Tephra2 codes in the directory ./src ,
		and RENAME the Tephra2 programme/executable to "Tephra2".
		(depending on the version of Tephra2 you are using, the Tephra2 executable might have different 
		names such as "tephra-2012"). 
		The existing files in ./src are Tephra2 installed on macOS. They might not be compatible to other operating systems.  	

	3. Put the observation mass/unit area data at some path as "observation.txt". The file should be a one 
		column vector (and no "," behind each value).
        	It should not have header/column name, and should be plain data.
 		Note: if thickness data is available, it needs to be transformed to
	     	mass/unit area.
 		See more details on Tephra2 websites. 
	4. Put the corresponding sample site coordinates and elevations as ./running_Tephra2/sites.csv.
		It should be a three column file without header/column names, and each element is separated by space (" ").
		It should include easting(UTM coordinate system), northing(UTM coordinate system), and elevation(m) of sample sites.  
		The sites' coordinates and elevations should be consistent with data in "observation.txt" (i.e., the first sample site
		corresponds to the first observed value in the two files).

***	5. Modify input files:
***		For both esp_input.csv and wind_input.csv, they have six columns.
			1st column: variable name (e.g., eruption source parameters and wind directions);
			2nd column: initial value of these variables to run the M-H algorithm;
			3rd column: whether a variable is known or not, and prior type of the variables to be estimated;
				It could have three values:
					"Gaussian"	the prior distribution for the variable of interest follows 
							a Gaussian distribution;
					"Uniform"	the prior distribution for the variable of interest follows 
							a Uniform distribution;
					"Fixed"		the value of this variable is known. Its value will be the 
							specified initial value (2nd column), and will not be changed
							during the implementation of the M-H algorithm.
			4th and 5th columns: parameters of the prior distribution:
				If the prior is "Gaussian", then it specifies the mean and standard deviation of the distribution;
				If the prior is "Uniform", then it species the min and max of the distribution.

			6th column: standard deviation/scale of the proposal distribution for each variable to be estimated.
***			*** For the variables with "Fixed" prior_type,
			*** their 4th-6th columns could be left blank.

			*** It should be noted that for this version of the algorithm, 
			*** we use four variables to specify the wind speed profile by assuming:
			*** wind direction does not change with elevation;
			*** wind speed increases from zero on the ground to a certain level linearly,
			*** and then starts to decrease with elevation until it reaches zero. 
			Therefore four variables are required to define the wind speed profile:
				(1) Wind direction;
				(2) Maximum wind speed;
				(3) Elevation that corresponds to the maximum wind speed;
				(4) Elevation where the wind speed reaches zero. 
***		Please note that currently users need to specify total eruption mass under log-scale (with the exponential constant as base).			


	5. Open spyder or other python IDE, load the script ./codes/run_1.7_codes.py, and follow the instructions within.
==========================================================================================================================================
	


