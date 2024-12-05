# fvDOMScatter
discrete ordinates with scattering (open beta)

This folder contains files and programs created under 
GNU GPL v3 License
by Mario Goddy Zuber 2024

(c) 2024 ETH Zurich, Mario Goddy Zuber

if you use any part of this work please cite the scientific contribution:
	
	M. Zuber, Dry Redox Reforming with Concentrated Solar Energy. 
	Thesis. ETH Zurich. Doctor of Sciences, 2024.

	

The model is implemented in OpenFOAM-v1906 (openfoam.com), which can be obtained from repository:
https://github.com/mgzuber-eth/fvDOMscatter
	
There is no dedicated user manual on this program. Please refer to the thesis 
of M. Zuber. In brief, to use the model, save the src files in a respective user
directory, and compile them via 'wmake'. To use the model create a new solver 
(e.g. chtMultiRegionFoam-custom) and include the following in the Make/options file:
-I(path_to_user_src_directory)/src/fvDOMScatter/lnInclude
-lfvDOMScatter

Simulation case files are included in the run directory. Run the simulations via './Allrun'.
The case files for a variation of the Heaslet and Warming problem with scattering phase function
are included. 
	
Note that the model is considered to be in a beta phase and requires further validation and work.
The majority of the scattering phase function/parameters are hard-coded into the src files.
	
