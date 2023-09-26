# swan-eps
Tools to facilitate ensemble forecasts using the SWAN model. Contains methods to download input data from ECMWF, 
write SWAN input files based on that data, and process the results. More specifically, this repository contains:

* **file_tools.py**:   tools to construct *.swn-input files, tools to construct shell scripts to run ensemble forecasts on Deltares computing clusters, and tools to download and manipulate data from ECMWF;
* **postprocess_tools.py**:  tools to plot SWAN output, convenient functions to manipulate strings, functions to predict computational costs of different types of ensemble forecasts, tools to interpolate SWAN output data in time and space and tools to estimate statistics from the results of SWAN ensemble forecasts;
* **construct_SWAN_case.ipynb** (not present yet):  example Jupyter notebook showing how one would go about constructing SWAN input based on ECMWF data using the tools in this repository;
* Example SWAN input and output to compare results:  wind input ensemble in boundary_conditions/meteo and boundary data in boundary_conditions/waves; bottom topography and output locations in geometry.

**FOLDER STRUCTURE**

In order for ensemble forecasts to work on your computer, your files have to be ordered in a very specific way. Inside a master folder, the following folders should be present:

* boundary_conditions:
    * meteo: contains ECMWF wind ensemble in netcdf-format compatible with SWAN; functions are available in file_tools.py to convert from ECMWF wind data to SWAN input;
    * waves: contains *.tpar files for nonstationary boundary condition data at all boundary points.
* geometry:
    * bottom topography file: *.bot-file defining the bottom level of your domain;
    * output_locations: contains files defining the output locations for which SWAN will return energy density spectrum data;
* EPSruns:
    * \{run1\} (Choose your own name)
        * input
            * *.swn-input files
            * shell script to run on a computing cluster
            * swaninit file (among other things, makes sure enough files can be loaded in)
        * output
    * \{run2\}
        * ...       
    * ...
 
**FILE-NAMING CONVENTIONS**

The code in this repository assumes certain conventions when it comes to the naming of files. These conventions are collected below.

* ***.tpar-files**: '\{chosen_name\}\_\{side\}\_\{index_starting_from_minimum\}.tpar'. Here, \{side\} is either 'w' (west), 'n' (north), 'e' (east) or 's' (south);
* **Wind ensemble members**: '\{chosen_name\}\{sample_index\}.nc'. Note that there is no underscore between \{chosen_name\} and \{sample_index\};
* **Standard ensemble forecast .swn-files**: '\{chosen_name\}\{sample_index\}.swn'. Again, there is no underscore between \{chosen_name\} and \{sample_index\};
* **Multilevel Monte Carlo ensemble forecast .swn-files**: '\{chosen_name\}\_level\{level_index\}\_sample\{sample_index\}.swn';
* **Multilevel multifidelity Monte Carlo ensemble forecast .swn-files**: '\{chosen_name\}\_level\{level_index\}_\{fidelity\}\_sample\{sample_index\}.swn'. Here, \{fidelity\} is either 'hi' or 'lo' to indicate high or low fidelity respectively;
* **Nested multilevel Monte Carlo ensemble forecast .swn-files**: '\{chosen_name\}\_olevel\{outer_level_index\}\_ilevel\{inner_level_index\}\_sample\{sample_index\}.swn'

**GENERAL NOTES**

* The code is written for an application to the North Sea domain. However, this can easily be changed in the code by changing the values of parameters in the Swanstring-class of 'file_tools.py'.
* To use these tools, you have to change the path variables at the start of the write_EPS_input-class in 'file_tools.py' to the paths that specify the location of your files. 

