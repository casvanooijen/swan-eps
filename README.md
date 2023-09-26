# swan-eps
Tools to facilitate ensemble forecasts using the SWAN model. Contains methods to download input data from ECMWF, 
write SWAN input files based on that data, and process the results. More specifically, this repository contains:

- **file_tools.py**:   tools to construct *.swn-input files, tools to construct shell scripts to run ensemble forecasts on Deltares computing clusters, and tools to download and manipulate data from ECMWF;
- **postprocess_tools.py**:  tools to plot SWAN output, convenient functions to manipulate strings, functions to predict computational costs of different types of ensemble forecasts, tools to interpolate SWAN output data in time and space and tools to estimate statistics from the results of SWAN ensemble forecasts;
- **construct_SWAN_case.ipynb** (not present yet):  example Jupyter notebook showing how one would go about constructing SWAN input based on ECMWF data using the tools in this repository;
- Example SWAN input and output to compare results:  wind input ensemble in boundary_conditions/meteo and boundary data in boundary_conditions/waves; bottom topography and output locations in geometry.

**Folder structure**
In order for ensemble forecasts to work on your computer, your files have to be ordered in a very specific way. A diagram is shown below:

                  boundary_conditions
                /
Master folder --- geometry
                \
                  EPSruns
