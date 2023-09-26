import xarray as xr
import numpy as np
import os
import postprocess_tools as ppt
from ecmwfapi import ECMWFDataServer


## MISCELLANEOUS ##

def write_swaninit(path, max_input_files=999, print_name="PRINT", input_name="INPUT"):
    """Writes a swaninit file to path.
    
    - path (str):               path to which swaninit will be written;
    - max_input_files (int):    maximum number of allowed input files in *.swn-files;
    - print_name (str):         name of the print-file that SWAN generates;
    - input_name (str):         name of the input-file that SWAN uses;
    """
    num_spaces = 5 - len(str(max_input_files))
    num_spaces_print = 40 - len(print_name)
    num_spaces_input = 40 - len(input_name)
    string = "    4                                   version of initialisation file\nDelft University of Technology          name of institute\n"
    string += f"    3                                   command file ref. number\n{input_name}" + num_spaces_input*" " + "input file name\n"
    string += f"    4                                   print file ref. number\n{print_name}" + num_spaces_print*" " + "print file name\n"
    string += f"    4                                   test file ref. number\n" + 40*" " + "test file name\n"
    string += f"    6                                   screen ref. number\n" + num_spaces * " " + f"{max_input_files}" + 35*" " + "highest file ref. number\n"
    string += f"$                                       comment identifier\n" + 43*" " + "TAB character\n"
    string += f"\\                                       dir sep char in input file\n/" + 39*" " + "dir sep char replacing previous one\n"
    string += f"    1                                   default time coding option"

    try:
        f = open(path + "swaninit", 'x')
        f.write(string)
    except FileExistsError:
        overwrite = input("Warning: overwriting other swaninit file. Continue [y/n]?")
        if overwrite = 'y':
            print("Overwriting")
            f = open(path + "swaninit", 'w')
            f.write(string)
        else:
            print("Cancelled")



## BUILD UP *.SWN INPUT FILES ##

class Swanstring(object):
    """This class is used to build up a SWAN input file for a North Sea application with *.swn-extension. 
    Use the methods input, boundary_conditions, physics, numerical_parameters and output to construct the text that the *.swn-file should 
    contain. Afterwards, use the method write_to_file to convert this text into a proper *.swn-file.
    
    - write_to_file:            writes the current string to a *.swn-file;
    - input:                    adds input paths and parameters to the text;
    - boundary_condition_tpar:  adds parametric, nonstationary boundary conditions to the text;
    - physics:                  adds physics commands to the text;
    - numerical_parameters:     defines the parameters for the numerical method and adds their definition to the text;
    - output:                   adds output commands and paths to the text.       
    """

    def __init__(self, fname, project, number, start_iso, end_iso, timestep):
        """Initialise Swanstring with non-stationary mode, spherical coordinates and nautical convention.
        
        - fname:        name of the SWAN input file (without *.swn-extension);
        - project:      name of the SWAN project;
        - number:       number of the run in this SWAN project;
        - start_iso:    start time in ISO-notation (e.g. October 10th 2001 is 20011010.000000);
        - end_iso:      end time in ISO-notation;
        - timestep:     amount of hours per timestep.        
        
        """
        self.fname = fname
        self.string = "!** pre-operational model version swan-ks (may 2023)\n!** for SWAN version 41.31A\n"+ f"!******************  HEADING  "+ \
            f"***************************\nPROJECT '{project}' '{number}'\n\n!******************  GENERAL  ***************************\n\nSET NAUTICAL\nSET LEVEL=0.000000;" \
                    + "\nSET CDCAP=0.002750\n\nMODE NONST\n\nCOORDINATES SPHERICAL CCM\n!******************  INPUT  ****************************\n"
        self.start_iso = start_iso
        self.end_iso = end_iso
        self.timestep = timestep # in hours
        
    
    def write_to_file(self, directory):
        try:
            f = open(directory + f"{self.fname}.swn", 'x')
            overwrite = True
        except FileExistsError:
            f = open(directory + f"{self.fname}.swn", 'w')
            overwrite = input('Warning: you will be overwriting an existing file. Continue [y/n]? ')
            if overwrite == 'y' or overwrite == 'Y' or overwrite == 'yes' or overwrite == 'Yes':
                overwrite = True
        if overwrite:
            f.write(self.string)


    def input(self, xnums, ynums, btmfile, watlfile=None, curfile=None, windfile=None, hotstart=None, spec_lowerbound=0.0300000):
        """Add INPUT commands to your text.
        
        - xnums (int):              number of grid cells in longitudinal direction;
        - ynums (int):              number of grid cells in latitudinal direction;
        - btmfile (str):            name of the bottom level file as one would write it in SWAN;
        - watlfile (str):           name of the water level file as one would write it in SWAN;
        - curfile (str):            name of the current (stroming) file as one would write it in SWAN;
        - windfile (str):           name of the wind input file as one would write it in SWAN;
        - hotstart (str):           name of the hotstart file as one would write it in SWAN;
        - spec_lowerbound (float):  lower bound of the spectral domain (default=0.03)

        """
        self.string += f"CGRID REGULAR -12 48 0.000000 21.000000 16.000000 {xnums} {ynums} CIRCLE 45 {spec_lowerbound} 0.600000\n\n"
        self.string += f"INP BOTTOM REGULAR  -12 48 0.000000 420 480 0.050000 0.033333 EXCE -999\n"
        self.string += f"READ BOTTOM 1. '{btmfile}' idla=5 FREE\n\n"

        if watlfile:
            self.string += f"INPGRID WLEV REGF '{watlfile}' NONSTAT {self.start_iso} {self.timestep} HR {self.end_iso}\n"
            self.string += f"READINP WLEV 1.0 '{watlfile}'\n\n"

        if curfile:
            self.string += f"INPGRID CUR REGF '{curfile}' NONSTAT {self.start_iso} {self.timestep} HR {self.end_iso}\n"
            self.string += f"READINP CUR 1.0 '{curfile}'\n\n"    

        if windfile:
            self.string += f"INPGRID WIND REGF '{windfile}' NONSTAT {self.start_iso} {self.timestep} HR {self.end_iso}\n"
            self.string += f"READINP WIND 1.0 '{windfile}'\n\n"

        if hotstart:
            self.string += f"INITIAL HOTSTART '{hotstart}' NETCDF\n\n"


    def boundary_condition_tpar(self, wavefield, tpar_fname, directory, tpar_swan_path, direction_param_type='PEAK'):
        """Adds (non-stationary) boundary condition commands to the text using *.tpar-files.

        - wavefield (xArray Dataset):   dataset from which the *.tpar-files were generated;
        - tpar_fname (str):             naming convention of the *.tpar-files (without '_{side}_{index}.tpar');
        - directory (str):              directory in which the *.tpar-files can be found;       
        - tpar_swan_path (str):         directory in which the *.tpar-files can be found as one would write it in a *.swn-input file;
        - direction_param_type (str):   type of direction parameter used ('PEAK' (default) or 'MEAN')        
        """
        
        self.string += "!******************  BOUNDARY CONDITIONS  **************\n\n"

        self.string += f"BOUND SHAPESPEC JONSWAP {direction_param_type} DSPR DEGREES\n\n" # ECMWF Dspr is given in degrees

        lon = wavefield['longitude']
        lat = wavefield['latitude']
        
        # West boundary
        # the following line assumes that at least one boundary point has boundary data
        self.string += f'BOUNDSPEC SIDE W CCW VARIABLE FILE &\n'
        for i in range(lat.values.shape[0]):
            dist = np.amax(lat.values) - lat.values[i]

            if os.path.isfile(directory+f'{tpar_fname}_w_{i}.tpar'):
                self.string += f"{dist} '{tpar_swan_path}{tpar_fname}_w_{i}.tpar' 1 & \n"

            if i == lat.values.shape[0] - 1:
                self.string = self.string[:-3] + '\n\n'
            
        
        # North boundary
        self.string += f'BOUNDSPEC SIDE N CCW VARIABLE FILE &\n'
        for i in range(lon.values.shape[0]):
            j = lon.values.shape[0] - 1 - i # make sure [len]-parameters appear in increasing order
            dist = np.amax(lon.values) - lon.values[j]

            if os.path.isfile(directory+f'{tpar_fname}_n_{j}.tpar'):
                self.string += f"{dist} '{tpar_swan_path}{tpar_fname}_n_{j}.tpar' 1 & \n"
                
            if i == lon.values.shape[0] - 1:
                self.string = self.string[:-3] + '\n\n'
            

        # East boundary
        self.string += f'BOUNDSPEC SIDE E CCW VARIABLE FILE &\n'
        for i in range(lat.values.shape[0]):
            j = lat.values.shape[0] - i - 1
            dist = lat.values[j] - np.amin(lat.values)

            if os.path.isfile(directory+f'{tpar_fname}_e_{j}.tpar'):
                self.string += f"{dist} '{tpar_swan_path}{tpar_fname}_e_{j}.tpar' 1 & \n"
            
            if i == lat.values.shape[0] - 1:
                self.string = self.string[:-3] + '\n\n'
    

        # South boundary
        self.string += f'BOUNDSPEC SIDE S CCW VARIABLE FILE &\n'
        for i in range(lon.values.shape[0]):
            dist = lon.values[i] - np.amin(lon.values)
            
            if os.path.isfile(directory+f'{tpar_fname}_s_{i}.tpar'):
                self.string += f"{dist} '{tpar_swan_path}{tpar_fname}_s_{i}.tpar' 1 & \n"
                
            if i == lon.values.shape[0] - 1:
                self.string = self.string[:-3] + '\n\n'


    def physics(self, generation=3):
        """Adds physics commands to the text. Parameter values are as is default in the North Sea model.
        
        - generation (int):     generation of SWAN (3 (default), 2, or 1).
        """
        self.string += "!******************  PHYSICS  **************************\n\n"
        if generation == 3:
            self.string += "GEN3 KOMEN\n"
        else:
            self.string += "GEN2\n"

        self.string += "WCAP KOMEN\nSSWELL\nQUAD    IQUAD=3\nFRIC    JONSWAP   CFJON=0.038\nBREA    CONST     ALPHA=1.0 GAMMA=0.73\n\n"


    def numerical_parameters(self, npnts=98, maxitns=20):
        """Adds commands pertaining to the numerical scheme this run of SWAN will use.
        
        - npnts (int):      percentage of grid points required to have converged (default=98);
        - maxitns (int):    maximum number of iterations (default=20).
        
        """

        self.string += "!******************  NUMERICAL PARAMETERS ****************\n\n"
        self.string += f"PROP BSBT\nNUMERIC ACCUR DREL=0.02 DHOVAL=0.02 DTOVAL=0.02 NPNTS={npnts} NONSTAT MXITNS={maxitns}\n\n"


    def output(self, output_locs_path, output_path, fname, output_vars = ['XP', "YP", "HSIG", "BOTLEV", "HSWELL", "TMM10", "TPS", "DIR", "DSPR", "WIND"]):
        """Add commands to decide what the output will look like (does not support *.tab-files). Output locations are as in 
        the default North Sea model. Output will be written to netcdf-format (*.nc).
        
        - output_locs_path (str):       path of output locations as one would write it in SWAN;
        - output_path (str):            path to which the output files should be written as one would write it in SWAN;
        - fname (str):                  naming convention of this SWAN run (without .nc);
        - output_vars (list):           list of output parameters one wishes to have SWAN compute. Default=['XP', 'YP',
                                        'HSIG', 'BOTLEV', 'HSWELL', 'TMM10', 'TPS', 'DIR', 'DSPR', 'WIND']
        
        """


        self.string += "!******************  OUTPUT ******************************\n"
        self.string += f"POINTS 'P1'    FILE '{output_locs_path}swan-ns_1.xyn'\nPOINTS 'P2'    FILE '{output_locs_path}swan-ns_2.xyn'\n" + \
            f"POINTS 'KS'    FILE '{output_locs_path}bndloc_kuststrook.xyn'\nPOINTS 'PM'    FILE '{output_locs_path}MDB.xyn'\n\n"
        
        self.string += f"BLOCK 'COMPGRID' NOHEAD  '{output_path}swan2d_{fname}.nc' LAYOUT 3 &\n "
        for var in output_vars:
            self.string += var + ' '
        self.string += f'OUT {self.start_iso} {self.timestep} HR\n\n'

        self.string += f"SPECOUT 'P1'    SPEC1D ABS '{output_path}swan-ns_spec1d_p1_{fname}.nc' 	    OUT {self.start_iso} {self.timestep} HR\n"
        self.string += f"SPECOUT 'P2'    SPEC2D ABS '{output_path}swan-ns_spec2d_p2_{fname}.nc'         OUT {self.start_iso} {self.timestep} HR\n"
        self.string += f"SPECOUT 'KS'    SPEC2D ABS '{output_path}swan-ns_spec2d_ks_{fname}.nc'         OUT {self.start_iso} {self.timestep} HR\n"
        self.string += f"SPECOUT 'PM'    SPEC2D ABS '{output_path}swan-ns_spec2d_pm_{fname}.nc'         OUT {self.start_iso} {self.timestep} HR\n\n"

        self.string += f"TEST 1 0\n\nCOMPUTE NONSTAT {self.start_iso} {self.timestep} HR {self.end_iso}\n\nSTOP"


## WRITE INPUT FILES FOR ENSEMBLE FORECASTS ##

class write_EPS_input:
    """This class contains methods to write input files for EPS runs with SWAN, using different methods. Your personal file paths may be defined in the beginning of this class.

    - era5wavefield (xArray dataset):                       ERA5 wave data from which boundary conditions are generated (one ensemble member or reanalysis);
    - bottom_path (str):                                    path to bottom file seen from the location of the *.swn input files, as one would write it in a *.swn-input file;
    - wind_path_without_filename_and_forecastdate (str):    partial path to wind ensemble files seen from the location of the *.swn input files, as one would write it in a *.swn-input file;
    - output_path (str):                                    path to location of eventual SWAN output files seen from the location of the *.swn input files, as one would write it in a *.swn-input file;
    - output_locs_path (str):                               path to files containing data of output locations seen from the location of the *.swn input files, as one would write it in a *.swn-input file;
    - tpar_filename (str):                                  name (without '_*_*.tpar') of *.tpar-boundary condition files;
    - tpar_swanpath (str):                                  path to *.tpar files seen from the location of the *.swn input files, as one would write it in a *.swn-input file;
    - tpar_directory (str):                                 full path to *.tpar files, which is compatible with file reading in python using 'open()'
    
    - standard_swn:     writes *.swn-input files for a standard Monte Carlo EPS run;
    - standard_sh:      writes a *.sh-file for a standard Monte Carlo EPS run on the Deltares h6 cluster;
    - ML_swn:           writes *.swn-input files for a multilevel Monte Carlo (ML) EPS run;
    - ML_sh:            writes a *.sh-file for a multilevel Monte Carlo EPS run on the Deltares h6 cluster;
    - MLMF_swn:         writes *.swn-input files for a multilevel multifidelity Monte Carlo(MLMF) EPS run;
    - MLMF_sh:          writes a *.sh-file for a multilevel multifidelity Monte Carlo EPS run on the Deltares h6 cluster;
    - NML_swn:          writes *.swn-input files for a nested multilevel Monte Carlo (NML) EPS run;
    - NML_sh:           writes a *.sh-file for a nested multilevel Monte Carlo EPS run on the Deltares h6 cluster.    
    """

    era5wavefield = xr.open_dataset("P:\\1230882-emodnet_hrsm\\vanCas\\Ensemble_Forecast\\boundary_conditions\\waves\\era5wavefield.nc", engine='netcdf4')
            
    # Define paths and filenames
    bottom_path = '../../geometry/swan-ns-j22_6-v1a_adjust.bot'
    wind_path_without_file_and_forecastdate = f'../../boundary_conditions/meteo/forecast_2013-12-' # forecast on December 4th 00:00:00
    output_path = '../output/'
    output_locs_path = '../../geometry/output_locations/'
    tpar_filename = 'stormxaver_boundary_wavedata'
    tpar_swanpath = '../../boundary_conditions/waves/boundary_condition_at_points/'
    tpar_directory = "P:\\1230882-emodnet_hrsm\\vanCas\\Ensemble_Forecast\\boundary_conditions\\waves\\boundary_condition_at_points\\"

    @staticmethod
    def standard_swn(path, fname, ensemble_size = 50):
        """Writes *.swn input files for a standard ensemble forecast, where ensemble members are ECMWF wind data.

        - path (str):           path to which the *.swn-files will be written;
        - fname (str):          name (without '_sample*.swn') of the *.swn-input files;
        - ensemble_size (int):  number of available ensemble members.
        """
        # Get input parameters interactively
        start_iso = input("Start time (ISO-notation): ")
        end_iso = input("End time (ISO-notation): ")
        timestep = float(input("Time step (def. 1): "))
        xnums = int(input("Number of cells in longitudinal direction (def. 420): "))
        ynums = int(input("Number of cells in latitudinal direction (def. 480): "))
        generation = int(input('Generation of SWAN (def. 3): '))
        spec_lowerbound = float(input("Lower bound of frequency domain (def. 0.03): "))
        maxitns = int(input("Maximum number of iterations per timestep (def. 20): "))

        for n in range(ensemble_size):
            swnstring = Swanstring(fname + f'_sample{n}', "SWAN-EPS", ppt.string.int_to_threedigitstring(n), start_iso, end_iso, timestep)
            swnstring.input(xnums, ynums, write_EPS_input.bottom_path, windfile=write_EPS_input.wind_path_without_file_and_forecastdate + f'{start_iso[6:8]}/stormxaver_wind_1204{n}.nc', spec_lowerbound=spec_lowerbound)
            swnstring.boundary_condition_tpar(write_EPS_input.era5wavefield, write_EPS_input.tpar_filename, write_EPS_input.tpar_directory, write_EPS_input.tpar_swanpath)
            #../../***** means that the shell script is two folder levels deeper than the folder boundary_conditions or geometry
            swnstring.physics(generation=generation)
            swnstring.numerical_parameters(maxitns=maxitns)
            swnstring.output(write_EPS_input.output_locs_path, write_EPS_input.output_path, fname + f'_sample{n}')
            swnstring.write_to_file(path)

    
    @staticmethod
    def standard_sh(path, fname, ensemble_size = 50):
        """Writes a *.sh-shell script to run a standard ensemble forecast on the Deltares h6-cluster. Path must be in the same folder as the *.swn input files and file must 
        have the same naming convention as the *.swn-files.
        
        - path (str):           path to which the *.sh-file will be written;
        - fname (str):  	    name of the *.sh-file without extension. Must be the same as the corresponding *.swn-files;
        - ensemble_size (int):  number of available ensemble members.
        """
        # Initial text
        filestring = f"#!/bin/sh\n\n#$ -cwd\n\n#$ -q normal-e3-c7\n\n#$ -N {fname}\n\n\n"
        # Load modules
        filestring += f"module load swan/41.31A.1_intel18.0.3_timers\nswan_omp_exe=swan_omp.exe\n\nexport OMP_NUM_THREADS=4"
        filestring += "echo ----------------------------------------------------------------------\n"
        filestring += "echo Run of\n"
        filestring += "echo $swan_omp_exe\n"
        filestring += "echo with OpenMP on linux-cluster.\n"
        filestring += "echo SGE_O_WORKDIR : $SGE_O_WORKDIR\n"
        filestring += "echo HOSTNAME : $HOSTNAME\n"
        filestring += "echo OMP_NUM_THREADS : $OMP_NUM_THREADS\n"
        filestring += "echo ----------------------------------------------------------------------\n\n"

        for n in range(ensemble_size):
            filestring += f"### Run {n}\n\n"
            filestring += f"cp {fname}_sample{n}.swn INPUT\n"
            filestring += "$swan_omp_exe\n"
            filestring += f"cp PRINT {fname}_sample{n}.prt\n\n"
            filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"

        f = open(path + fname + '.sh', 'x')
        f.write(filestring)


    @staticmethod
    def ML_swn(path, fname, N_list):
        """Writes *.swn-input files for an ensemble forecast using the multilevel Monte Carlo method. Different approximation types are available
        
        - path (str):       path to which the *.swn-files will be written;
        - fname (str):      name of the *.swn-files without '_level*_sample*.swn';
        - N_list (list):    values of N used in the multilevel Monte Carlo method.
        """
        # Select approximation type
        valid = False
        while not valid:
            approx_type = input("Select approximation type:\nspatial_resolution3 (three different spatial resolutions)\n")
            if approx_type in ['spatial_resolution3']:
                valid = True

        if approx_type == 'spatial_resolution3':
            # Get input parameters interactively
            start_iso = input("Start time (ISO-notation) (def. 20131204.000000): ")
            end_iso = input("End time (ISO-notation) (def. 20131209.000000): ")
            timestep = float(input("Time step (def. 1): "))
            generation = int(input('Generation of SWAN (def. 3): '))
            spec_lowerbound = float(input("Lower bound of frequency domain (def. 0.03): "))
            maxitns = int(input("Maximum number of iterations per timestep (def. 20): "))
            # Skip number of grid cells, because these are predefined in this type of multilevel approximation

            # Open ERA5 wave data; necessary to write boundary conditions;
            era5wavefield = write_EPS_input.era5wavefield
            
            # Define paths and filenames
            bottom_path = write_EPS_input.bottom_path
            wind_path_without_file = write_EPS_input.wind_path_without_file_and_forecastdate + f'{start_iso[6:8]}/'
            output_path = write_EPS_input.output_path
            output_locs_path = write_EPS_input.output_locs_path
            tpar_filename = write_EPS_input.tpar_filename
            tpar_swanpath = write_EPS_input.tpar_swanpath
            tpar_directory = write_EPS_input.tpar_directory

            # LEVEL 0
            for n in range(N_list[0]):
                sample_number = n
                swnstring = Swanstring(fname + f'_level0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(n), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file + f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=generation)
                swnstring.numerical_parameters(maxitns=maxitns)
                swnstring.output(output_locs_path, output_path, fname + f'_level0_sample{sample_number}')
                swnstring.write_to_file(path)

            # HIGHER LEVEL CORRECTIONS
            for l in range(1, 3):
                for n in range(N_list[l]):
                    sample_number = sum(N_list[:l]) + n
                    run_number = N_list[0] + 2*sum(N_list[1:l]) + n

                    swnstring_lo = Swanstring(fname+f'_level{l-1}_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_number), start_iso, end_iso, timestep)
                    swnstring_hi = Swanstring(fname+f'_level{l}_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_number + N_list[l]), start_iso, end_iso, timestep)

                    swnstring_lo.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_hi.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)

                    swnstring_lo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_hi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)

                    swnstring_lo.physics(generation=generation)
                    swnstring_hi.physics(generation=generation)

                    swnstring_lo.numerical_parameters(maxitns=maxitns)
                    swnstring_hi.numerical_parameters(maxitns=maxitns)
                    
                    swnstring_lo.output(output_locs_path, output_path, fname + f'level{l-1}_sample{sample_number}')
                    swnstring_hi.output(output_locs_path, output_path, fname + f'level{l}_sample{sample_number}')

                    swnstring_lo.write_to_file(path)
                    swnstring_hi.write_to_file(path)        


    @staticmethod
    def ML_sh(path, fname, N_list):
        """Writes a *.sh-shell script to run a multilevel Monte Carlo ensemble forecast on the Deltares h6 cluster. Path must be in the same folder as the *.swn input files and file must 
        have the same naming convention as the *.swn-files.
        
        - path (str):       path to which the *.sh-file will be written;
        - fname (str):      name of the *.sh-file without extension. Must be the same as the corresponding *.swn-files;
        - N_list (str):     values of N used in the multilevel Monte Carlo method.
        """
        # Initial text
        filestring = f"#!/bin/sh\n\n#$ -cwd\n\n#$ -q normal-e3-c7\n\n#$ -N {fname}\n\n\n"
        # Load modules
        filestring += f"module load swan/41.31A.1_intel18.0.3_timers\nswan_omp_exe=swan_omp.exe\n\nexport OMP_NUM_THREADS=4"
        filestring += "echo ----------------------------------------------------------------------\n"
        filestring += "echo Run of\n"
        filestring += "echo $swan_omp_exe\n"
        filestring += "echo with OpenMP on linux-cluster.\n"
        filestring += "echo SGE_O_WORKDIR : $SGE_O_WORKDIR\n"
        filestring += "echo HOSTNAME : $HOSTNAME\n"
        filestring += "echo OMP_NUM_THREADS : $OMP_NUM_THREADS\n"
        filestring += "echo ----------------------------------------------------------------------\n\n"

        # LEVEL 0
        for n in range(N_list[0]):
            sample_number = n
            filestring += f"### LEVEL 0, SAMPLE {sample_number}\n\n"
            filestring += f"cp {fname}_level0_sample{sample_number}.swn INPUT\n"
            filestring += "$swan_omp_exe\n"
            filestring += f"cp PRINT {fname}_level0_sample{sample_number}.prt\n"
            filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"

        # HIGHER LEVEL CORRECTIONS
        for l in range(1, len(N_list)):
            for n in range(N_list[l]):
                sample_number = sum(N_list[:l]) + n
                filestring += f"### LEVEL {l-1}, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_level{l-1}_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_level{l-1}_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
            for n in range(N_list[l]):
                sample_number = sum(N_list[:l]) + n
                filestring += f"### LEVEL {l}, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_level{l}_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_level{l}_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
        
        f = open(path + fname + '.sh', 'x')
        f.write(filestring)


    @staticmethod
    def MLMF_swn(path, fname, N_list, r_list):
        """Writes *.swn-input files for an ensemble forecast using the multilevel multifidelity Monte Carlo method.
        
        - path (str):       path to which the *.swn-files will be written;
        - fname (str):      name of the *.swn-files without '_level*_*_sample*.swn';
        - N_list (list):    values of N^(HF) used in the multilevel multifidelity Monte Carlo method;
        - r_list (list):    values of r used in the multilevel multifidelity Monte Carlo method.
        """
        # Select approximation type
        valid = False
        while not valid:
            approx_type = input("You are trying to construct *.swn-input files for a multilevel multifidelity Monte Carlo EPS simulation. " 
                                +"Select the approximation type:\n lspace_ffreqitns (levels: spatial resolution, fidelity: smaller frequency domain and more tolerant stopping criterion (maximum number of iterations is lower));\n"
                                +"lspace_fgen2 (levels: spatial resolution, fidelity: swan generation 2);\n")
            if approx_type in ['lspace_ffreqitns', 'lspace_fgen2']:
                valid = True
        
        if approx_type == "lspace_ffreqitns":
            # Get input parameters interactively
            start_iso = input("Start time (ISO-notation) (def. 20131204.000000): ")
            end_iso = input("End time (ISO-notation) (def. 20131209.000000): ")
            timestep = float(input("Time step (def. 1): "))
            generation = int(input('Generation of SWAN (def. 3): '))
        	# Skip number of grid cells, maximum iteration and frequency domain, because these are fixed by the approximation method

            # Open ERA5 wave data; necessary to write boundary conditions;
            era5wavefield = write_EPS_input.era5wavefield
            
            # Define paths and filenames
            bottom_path = write_EPS_input.bottom_path
            wind_path_without_file = write_EPS_input.wind_path_without_file_and_forecastdate + f'{start_iso[6:8]}/'
            output_path = write_EPS_input.output_path
            output_locs_path = write_EPS_input.output_locs_path
            tpar_filename = write_EPS_input.tpar_filename
            tpar_swanpath = write_EPS_input.tpar_swanpath
            tpar_directory = write_EPS_input.tpar_directory

            # Level 0
            # high-fidelity
            run_counter = 1

            for n in range(N_list[0]):
                sample_number = n
                swnstring = Swanstring(fname + f'_level0_hi_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.03)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=generation)
                swnstring.numerical_parameters(maxitns=20)
                swnstring.output(output_locs_path, output_path, fname + f'_level0_hi_sample{sample_number}')
                swnstring.write_to_file(path)
                run_counter += 1
            
            # low-fidelity
            for sample_number in ppt.estimate.used_sample_indices_MLMF(N_list, r_list, 0):
                swnstring = Swanstring(fname + f'_level0_lo_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=generation)
                swnstring.numerical_parameters(maxitns=3)
                swnstring.output(output_locs_path, output_path, fname + f'_level0_lo_sample{sample_number}')
                swnstring.write_to_file(path)
                run_counter += 1

            # Higher level corrections
            for l in range(1, 3):
                # high-fidelity
                for n in range(N_list[l]):
                    sample_number = sum(N_list[:l]) + n

                    swnstring_lowlevel = Swanstring(fname+f'_level{l-1}_hi_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_highlevel = Swanstring(fname+f'_level{l}_hi_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter + 1), start_iso, end_iso, timestep)

                    swnstring_lowlevel.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.03)
                    swnstring_highlevel.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.03)

                    swnstring_lowlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_highlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)

                    swnstring_lowlevel.physics(generation=generation)
                    swnstring_highlevel.physics(generation=generation)

                    swnstring_lowlevel.numerical_parameters(maxitns=20)
                    swnstring_highlevel.numerical_parameters(maxitns=20)
                    
                    swnstring_lowlevel.output(output_locs_path, output_path, fname + f'level{l-1}_hi_sample{sample_number}')
                    swnstring_highlevel.output(output_locs_path, output_path, fname + f'level{l}_hi_sample{sample_number}')

                    swnstring_lowlevel.write_to_file(path)
                    swnstring_highlevel.write_to_file(path)

                    run_counter += 2

                # low-fidelity
                for sample_number in ppt.estimate.used_sample_indices_MLMF(N_list, r_list, l):
                    swnstring = Swanstring(fname + f'_level{l}_lo_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring.input(105*(2**l), 120*(2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                    swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring.physics(generation=generation)
                    swnstring.numerical_parameters(maxitns=3)
                    swnstring.output(output_locs_path, output_path, fname + f'_level{l}_lo_sample{sample_number}')
                    swnstring.write_to_file(path)
                    run_counter += 1
                    
        if approx_type == "lspace_fgen2":
            # Get input parameters interactively
            start_iso = input("Start time (ISO-notation) (def. 20131204.000000): ")
            end_iso = input("End time (ISO-notation) (def. 20131209.000000): ")
            timestep = float(input("Time step (def. 1): "))
            maxitns = int(input('Maximum number of iterations per timestep (def. 20): '))
            spec_lowerbound = float(input('Lower bound of frequency domain (def. 0.03): '))
        	# Skip number of grid cells, and generation, because these are fixed by the approximation method

            # Open ERA5 wave data; necessary to write boundary conditions;
            era5wavefield = write_EPS_input.era5wavefield
            
            # Define paths and filenames
            bottom_path = write_EPS_input.bottom_path
            wind_path_without_file = write_EPS_input.wind_path_without_file_and_forecastdate + f'{start_iso[6:8]}/'
            output_path = write_EPS_input.output_path
            output_locs_path = write_EPS_input.output_locs_path
            tpar_filename = write_EPS_input.tpar_filename
            tpar_swanpath = write_EPS_input.tpar_swanpath
            tpar_directory = write_EPS_input.tpar_directory

            # Level 0
            # high-fidelity
            run_counter = 1

            for n in range(N_list[0]):
                sample_number = n
                swnstring = Swanstring(fname + f'_level0_hi_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=3)
                swnstring.numerical_parameters(maxitns=maxitns)
                swnstring.output(output_locs_path, output_path, fname + f'_level0_hi_sample{sample_number}')
                swnstring.write_to_file(path)
                run_counter += 1
            
            # low-fidelity
            for sample_number in ppt.estimate.used_sample_indices_MLMF(N_list, r_list, 0):
                swnstring = Swanstring(fname + f'_level0_lo_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=2)
                swnstring.numerical_parameters(maxitns=maxitns)
                swnstring.output(output_locs_path, output_path, fname + f'_level0_lo_sample{sample_number}')
                swnstring.write_to_file(path)
                run_counter += 1

            # Higher level corrections
            for l in range(1, 3):
                # high-fidelity
                for n in range(N_list[l]):
                    sample_number = sum(N_list[:l]) + n

                    swnstring_lowlevel = Swanstring(fname+f'_level{l-1}_hi_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_highlevel = Swanstring(fname+f'_level{l}_hi_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter + 1), start_iso, end_iso, timestep)

                    swnstring_lowlevel.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_highlevel.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)

                    swnstring_lowlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_highlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)

                    swnstring_lowlevel.physics(generation=3)
                    swnstring_highlevel.physics(generation=3)

                    swnstring_lowlevel.numerical_parameters(maxitns=maxitns)
                    swnstring_highlevel.numerical_parameters(maxitns=maxitns)
                    
                    swnstring_lowlevel.output(output_locs_path, output_path, fname + f'level{l-1}_hi_sample{sample_number}')
                    swnstring_highlevel.output(output_locs_path, output_path, fname + f'level{l}_hi_sample{sample_number}')

                    swnstring_lowlevel.write_to_file(path)
                    swnstring_highlevel.write_to_file(path)

                    run_counter += 2

                # low-fidelity
                for sample_number in ppt.estimate.used_sample_indices_MLMF(N_list, r_list, l):
                    swnstring = Swanstring(fname + f'_level{l}_lo_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring.input(105*(2**l), 120*(2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring.physics(generation=2)
                    swnstring.numerical_parameters(maxitns=maxitns)
                    swnstring.output(output_locs_path, output_path, fname + f'_level{l}_lo_sample{sample_number}')
                    swnstring.write_to_file(path)
                    run_counter += 1


    @staticmethod
    def MLMF_sh(path, fname, N_list, r_list):
        """Writes a *.sh-shell script to run a multilevel multifidelity Monte Carlo ensemble forecast on the Deltares h6 cluster. Path must be in the same folder as the *.swn input files and file
        must have the same naming convention as the *.swn-files.
        
        - path (str):       path to which the *.sh-file will be written;
        - fname (str):      name of the *.sh-file without extension. Must be the same as the corresponding *.swn-files;
        - N_list (list):    values of N^(HF) used in the multilevel multifidelity Monte Carlo method;
        - r_list (list):    values of r used in the multilevel multifidelity Monte Carlo method.
        """
        # Initial text
        filestring = f"#!/bin/sh\n\n#$ -cwd\n\n#$ -q normal-e3-c7\n\n#$ -N {fname}\n\n\n"
        # Load modules
        filestring += f"module load swan/41.31A.1_intel18.0.3_timers\nswan_omp_exe=swan_omp.exe\n\nexport OMP_NUM_THREADS=4"
        filestring += "echo ----------------------------------------------------------------------\n"
        filestring += "echo Run of\n"
        filestring += "echo $swan_omp_exe\n"
        filestring += "echo with OpenMP on linux-cluster.\n"
        filestring += "echo SGE_O_WORKDIR : $SGE_O_WORKDIR\n"
        filestring += "echo HOSTNAME : $HOSTNAME\n"
        filestring += "echo OMP_NUM_THREADS : $OMP_NUM_THREADS\n"
        filestring += "echo ----------------------------------------------------------------------\n\n"

        for n in range(N_list[0]):
            sample_number = n
            filestring += f"### LEVEL 0, HIGH-FIDELITY, SAMPLE {sample_number}\n\n"
            filestring += f"cp {fname}_level0_hi_sample{sample_number}.swn INPUT\n"
            filestring += "$swan_omp_exe\n"
            filestring += f"cp PRINT {fname}_level0_hi_sample{sample_number}.prt\n"
            filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
        for sample_number in ppt.estimate.used_sample_indices_MLMF(N_list, r_list, 0):
            filestring += f"### LEVEL 0, LOW-FIDELITY, SAMPLE {sample_number}\n\n"
            filestring += f"cp {fname}_level0_lo_sample{sample_number}.swn INPUT\n"
            filestring += "$swan_omp_exe\n"
            filestring += f"cp PRINT {fname}_level0_lo_sample{sample_number}.prt\n"
            filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"

        for l in range(1, len(N_list)):
            for n in range(N_list[l]):
                sample_number = sum(N_list[:l]) + n
                filestring += f"### LEVEL {l-1}, HIGH-FIDELITY, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_level{l-1}_hi_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_level{l-1}_hi_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
            for n in range(N_list[l]):
                sample_number = sum(N_list[:l]) + n
                filestring += f"### LEVEL {l}, HIGH-FIDELITY, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_level{l}_hi_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_level{l}_hi_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
            for sample_number in ppt.estimate.used_sample_indices_MLMF(N_list, r_list, l):
                filestring += f"### LEVEL {l}, LOW-FIDELITY, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_level{l}_lo_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_level{l}_lo_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"

        f = open(path + fname + '.sh', 'x')
        f.write(filestring)


    @staticmethod
    def NML_swn(path, fname, M_array):
        """Writes *.swn-input files for ensemble forecasts using the nested multilevel Monte Carlo method.
        
        - path (str):           path to which the *.swn-files will be written;
        - fname (str):          name of the *.swn-files without '_olevel*_ilevel*_sample*.swn';
        - M_array (np.ndarray): values of M used in the nested multilevel Monte Carlo method.
        """
        # Select approximation type
        valid = False
        while not valid:
            approx_type = input("You are trying to construct *.swn-input files for a nested multilevel Monte Carlo EPS simulation." + 
                                "Select the approximation type.\n" + 
                                "ospace_ifreqitns (outer levels: spatial resolution, inner level: smaller frequency domain and more tolerant stopping criterion);\n" + 
                                "ospace_igen2 (outer levels: spatial resolution, inner level: swan generation 2);\n")
            if approx_type in ['ospace_ifreqitns', 'ospace_igen2']:
                valid = True

        if approx_type == 'ospace_ifreqitns':
            # Get input parameters interactively
            start_iso = input("Start time (ISO-notation) (def. 20131204.000000): ")
            end_iso = input("End time (ISO-notation) (def. 20131209.000000): ")
            timestep = float(input("Time step (def. 1): "))
            generation = int(input('Generation of SWAN (def. 3): '))
        	# Skip number of grid cells, maximum iteration and frequency domain, because these are fixed by the approximation method

            # Open ERA5 wave data; necessary to write boundary conditions;
            era5wavefield = write_EPS_input.era5wavefield
            
            # Define paths and filenames
            bottom_path = write_EPS_input.bottom_path
            wind_path_without_file = write_EPS_input.wind_path_without_file_and_forecastdate + f'{start_iso[6:8]}/'
            output_path = write_EPS_input.output_path
            output_locs_path = write_EPS_input.output_locs_path
            tpar_filename = write_EPS_input.tpar_filename
            tpar_swanpath = write_EPS_input.tpar_swanpath
            tpar_directory = write_EPS_input.tpar_directory

            run_counter = 1

            # Outer level 0
            for n in range(M_array[0, 0]):
                sample_number = n
                swnstring = Swanstring(fname + f'_olevel0_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=generation)
                swnstring.numerical_parameters(maxitns=3)
                swnstring.output(output_locs_path, output_path, fname + f'_olevel0_ilevel0_sample{sample_number}')
                swnstring.write_to_file(path)
                run_counter += 1
            
            for n in range(M_array[0, 1]):
                sample_number = np.sum(M_array[0, :1]) + n
                swnstring_lo = Swanstring(fname + f'_olevel0_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring_lo.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                swnstring_lo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring_lo.physics(generation=generation)
                swnstring_lo.numerical_parameters(maxitns=3)
                swnstring_lo.output(output_locs_path, output_path, fname + f'_olevel0_ilevel0_sample{sample_number}')
                swnstring_lo.write_to_file(path)
                run_counter += 1

                swnstring_hi = Swanstring(fname + f'_olevel0_ilevel1_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring_hi.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.03)
                swnstring_hi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring_hi.physics(generation=generation)
                swnstring_hi.numerical_parameters(maxitns=20)
                swnstring_hi.output(output_locs_path, output_path, fname + f'_olevel0_ilevel1_sample{sample_number}')
                swnstring_hi.write_to_file(path)
                run_counter += 1

            # Higher outer level corrections
            for l in range(1, 3):
                for n in range(M_array[l, 0]):
                    sample_number = np.sum(M_array[:l, :]) + n

                    swnstring_lowlevel = Swanstring(fname + f'_olevel{l-1}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_lowlevel.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                    swnstring_lowlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_lowlevel.physics(generation=generation)
                    swnstring_lowlevel.numerical_parameters(maxitns=3)
                    swnstring_lowlevel.output(output_locs_path, output_path, fname + f'_olevel{l-1}_ilevel0_sample{sample_number}')
                    swnstring_lowlevel.write_to_file(path)
                    run_counter += 1

                    swnstring_highlevel = Swanstring(fname + f'_olevel{l}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_highlevel.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                    swnstring_highlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_highlevel.physics(generation=generation)
                    swnstring_highlevel.numerical_parameters(maxitns=3)
                    swnstring_highlevel.output(output_locs_path, output_path, fname + f'_olevel{l}_ilevel0_sample{sample_number}')
                    swnstring_highlevel.write_to_file(path)
                    run_counter += 1

                for n in range(M_array[l, 1]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :1]) + n

                    swnstring_olo_ilo = Swanstring(fname + f'_olevel{l-1}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_olo_ilo.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                    swnstring_olo_ilo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_olo_ilo.physics(generation=generation)
                    swnstring_olo_ilo.numerical_parameters(maxitns=3)
                    swnstring_olo_ilo.output(output_locs_path, output_path, fname + f'_olevel{l-1}_ilevel0_sample{sample_number}')
                    swnstring_olo_ilo.write_to_file(path)
                    run_counter += 1

                    swnstring_ohi_ilo = Swanstring(fname + f'_olevel{l}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_ohi_ilo.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.06)
                    swnstring_ohi_ilo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_ohi_ilo.physics(generation=generation)
                    swnstring_ohi_ilo.numerical_parameters(maxitns=3)
                    swnstring_ohi_ilo.output(output_locs_path, output_path, fname + f'_olevel{l}_ilevel0_sample{sample_number}')
                    swnstring_ohi_ilo.write_to_file(path)
                    run_counter += 1

                    swnstring_olo_ihi = Swanstring(fname + f'_olevel{l-1}_ilevel1_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_olo_ihi.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.03)
                    swnstring_olo_ihi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_olo_ihi.physics(generation=generation)
                    swnstring_olo_ihi.numerical_parameters(maxitns=20)
                    swnstring_olo_ihi.output(output_locs_path, output_path, fname + f'_olevel{l-1}_ilevel1_sample{sample_number}')
                    swnstring_olo_ihi.write_to_file(path)
                    run_counter += 1

                    swnstring_ohi_ihi = Swanstring(fname + f'_olevel{l}_ilevel1_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_ohi_ihi.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=0.03)
                    swnstring_ohi_ihi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_ohi_ihi.physics(generation=generation)
                    swnstring_ohi_ihi.numerical_parameters(maxitns=20)
                    swnstring_ohi_ihi.output(output_locs_path, output_path, fname + f'_olevel{l}_ilevel1_sample{sample_number}')
                    swnstring_ohi_ihi.write_to_file(path)
                    run_counter += 1

        if approx_type == 'ospace_igen2':
            # Get input parameters interactively
            start_iso = input("Start time (ISO-notation) (def. 20131204.000000): ")
            end_iso = input("End time (ISO-notation) (def. 20131209.000000): ")
            timestep = float(input("Time step (def. 1): "))
            maxitns = int(input('Maximum number of iterations per timestep (def. 20): '))
            spec_lowerbound = float(input('Lower bound of frequency domain (def. 0.03): '))
        	# Skip number of grid cells, maximum iteration and frequency domain, because these are fixed by the approximation method

            # Open ERA5 wave data; necessary to write boundary conditions;
            era5wavefield = write_EPS_input.era5wavefield
            
            # Define paths and filenames
            bottom_path = write_EPS_input.bottom_path
            wind_path_without_file = write_EPS_input.wind_path_without_file_and_forecastdate + f'{start_iso[6:8]}/'
            output_path = write_EPS_input.output_path
            output_locs_path = write_EPS_input.output_locs_path
            tpar_filename = write_EPS_input.tpar_filename
            tpar_swanpath = write_EPS_input.tpar_swanpath
            tpar_directory = write_EPS_input.tpar_directory

            run_counter = 1

            # Outer level 0
            for n in range(M_array[0, 0]):
                sample_number = n
                swnstring = Swanstring(fname + f'_olevel0_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                swnstring.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring.physics(generation=2)
                swnstring.numerical_parameters(maxitns=maxitns)
                swnstring.output(output_locs_path, output_path, fname + f'_olevel0_ilevel0_sample{sample_number}')
                swnstring.write_to_file(path)
                run_counter += 1
            
            for n in range(M_array[0, 1]):
                sample_number = np.sum(M_array[0, :1]) + n
                swnstring_lo = Swanstring(fname + f'_olevel0_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring_lo.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                swnstring_lo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring_lo.physics(generation=2)
                swnstring_lo.numerical_parameters(maxitns=maxitns)
                swnstring_lo.output(output_locs_path, output_path, fname + f'_olevel0_ilevel0_sample{sample_number}')
                swnstring_lo.write_to_file(path)
                run_counter += 1

                swnstring_hi = Swanstring(fname + f'_olevel0_ilevel1_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                swnstring_hi.input(105, 120, bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                swnstring_hi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                swnstring_hi.physics(generation=3)
                swnstring_hi.numerical_parameters(maxitns=maxitns)
                swnstring_hi.output(output_locs_path, output_path, fname + f'_olevel0_ilevel1_sample{sample_number}')
                swnstring_hi.write_to_file(path)
                run_counter += 1

            # Higher outer level corrections
            for l in range(1, 3):
                for n in range(M_array[l, 0]):
                    sample_number = np.sum(M_array[:l, :]) + n

                    swnstring_lowlevel = Swanstring(fname + f'_olevel{l-1}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_lowlevel.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_lowlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_lowlevel.physics(generation=2)
                    swnstring_lowlevel.numerical_parameters(maxitns=maxitns)
                    swnstring_lowlevel.output(output_locs_path, output_path, fname + f'_olevel{l-1}_ilevel0_sample{sample_number}')
                    swnstring_lowlevel.write_to_file(path)
                    run_counter += 1

                    swnstring_highlevel = Swanstring(fname + f'_olevel{l}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_highlevel.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_highlevel.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_highlevel.physics(generation=2)
                    swnstring_highlevel.numerical_parameters(maxitns=maxitns)
                    swnstring_highlevel.output(output_locs_path, output_path, fname + f'_olevel{l}_ilevel0_sample{sample_number}')
                    swnstring_highlevel.write_to_file(path)
                    run_counter += 1

                for n in range(M_array[l, 1]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :1]) + n

                    swnstring_olo_ilo = Swanstring(fname + f'_olevel{l-1}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_olo_ilo.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_olo_ilo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_olo_ilo.physics(generation=2)
                    swnstring_olo_ilo.numerical_parameters(maxitns=maxitns)
                    swnstring_olo_ilo.output(output_locs_path, output_path, fname + f'_olevel{l-1}_ilevel0_sample{sample_number}')
                    swnstring_olo_ilo.write_to_file(path)
                    run_counter += 1

                    swnstring_ohi_ilo = Swanstring(fname + f'_olevel{l}_ilevel0_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_ohi_ilo.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_ohi_ilo.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_ohi_ilo.physics(generation=2)
                    swnstring_ohi_ilo.numerical_parameters(maxitns=maxitns)
                    swnstring_ohi_ilo.output(output_locs_path, output_path, fname + f'_olevel{l}_ilevel0_sample{sample_number}')
                    swnstring_ohi_ilo.write_to_file(path)
                    run_counter += 1

                    swnstring_olo_ihi = Swanstring(fname + f'_olevel{l-1}_ilevel1_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_olo_ihi.input(105 * (2**(l-1)), 120 * (2**(l-1)), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_olo_ihi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_olo_ihi.physics(generation=3)
                    swnstring_olo_ihi.numerical_parameters(maxitns=maxitns)
                    swnstring_olo_ihi.output(output_locs_path, output_path, fname + f'_olevel{l-1}_ilevel1_sample{sample_number}')
                    swnstring_olo_ihi.write_to_file(path)
                    run_counter += 1

                    swnstring_ohi_ihi = Swanstring(fname + f'_olevel{l}_ilevel1_sample{sample_number}', 'SWAN-EPS', ppt.string.int_to_threedigitstring(run_counter), start_iso, end_iso, timestep)
                    swnstring_ohi_ihi.input(105 * (2**l), 120 * (2**l), bottom_path, windfile=wind_path_without_file+f'stormxaver_wind_1204{sample_number}.nc', spec_lowerbound=spec_lowerbound)
                    swnstring_ohi_ihi.boundary_condition_tpar(era5wavefield, tpar_filename, tpar_directory, tpar_swanpath)
                    swnstring_ohi_ihi.physics(generation=3)
                    swnstring_ohi_ihi.numerical_parameters(maxitns=maxitns)
                    swnstring_ohi_ihi.output(output_locs_path, output_path, fname + f'_olevel{l}_ilevel1_sample{sample_number}')
                    swnstring_ohi_ihi.write_to_file(path)
                    run_counter += 1


    @staticmethod
    def NML_sh(path, fname, M_array):
        """Writes a *.sh-shell script to run a nested multilevel Monte Carlo ensemble forecast. Path must be in the same folder as the *.swn input files and the file
        must have the same naming convention as the *.swn-files.
        
        - path (str):           path to which the *.sh-file will be written;
        - fname (str):          name of the *.sh-file without extension. Must be the same as the corresponding *.swn-files;
        - M_array (np.ndarray): values of M used in the nested multilevel Monte Carlo method.
        """
        num_outer_levels = M_array.shape[0]
        num_inner_levels = M_array.shape[1]

        # Initial text
        filestring = f"#!/bin/sh\n\n#$ -cwd\n\n#$ -q normal-e3-c7\n\n#$ -N {fname}\n\n\n"
        # Load modules
        filestring += f"module load swan/41.31A.1_intel18.0.3_timers\nswan_omp_exe=swan_omp.exe\n\nexport OMP_NUM_THREADS=4"
        filestring += "echo ----------------------------------------------------------------------\n"
        filestring += "echo Run of\n"
        filestring += "echo $swan_omp_exe\n"
        filestring += "echo with OpenMP on linux-cluster.\n"
        filestring += "echo SGE_O_WORKDIR : $SGE_O_WORKDIR\n"
        filestring += "echo HOSTNAME : $HOSTNAME\n"
        filestring += "echo OMP_NUM_THREADS : $OMP_NUM_THREADS\n"
        filestring += "echo ----------------------------------------------------------------------\n\n"

        # Outer level 0
        for n in range(M_array[0,0]):
            sample_number = n
            filestring += f"### OUTER LEVEL 0, INNER LEVEL 0, SAMPLE {sample_number}\n\n"
            filestring += f"cp {fname}_olevel0_ilevel0_sample{sample_number}.swn INPUT\n"
            filestring += "$swan_omp_exe\n"
            filestring += f"cp PRINT {fname}_olevel0_ilevel0_sample{sample_number}.prt\n"
            filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"

        for k in range(1, num_inner_levels):
            for n in range(M_array[0, k]):
                sample_number = np.sum(M_array[0, :k]) + n
                filestring += f"### OUTER LEVEL 0, INNER LEVEL {k-1}, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_olevel0_ilevel{k-1}_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_olevel0_ilevel{k-1}_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
            for n in range(M_array[0, k]):
                sample_number = np.sum(M_array[0, :k]) + n
                filestring += f"### OUTER LEVEL 0, INNER LEVEL {k}, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_olevel0_ilevel{k}_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_olevel0_ilevel{k}_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"

        # Higher outer level corrections
        for l in range(1, num_outer_levels):
            for n in range(M_array[l, 0]):
                sample_number = np.sum(M_array[:l, :]) + n
                filestring += f"### OUTER LEVEL {l-1}, INNER LEVEL 0, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_olevel{l-1}_ilevel0_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_olevel{l-1}_ilevel0_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
            for n in range(M_array[l, 0]):
                sample_number = np.sum(M_array[:l, :]) + n
                filestring += f"### OUTER LEVEL {l}, INNER LEVEL 0, SAMPLE {sample_number}\n\n"
                filestring += f"cp {fname}_olevel{l}_ilevel0_sample{sample_number}.swn INPUT\n"
                filestring += "$swan_omp_exe\n"
                filestring += f"cp PRINT {fname}_olevel{l}_ilevel0_sample{sample_number}.prt\n"
                filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
            
            for k in range(1, num_inner_levels):
                for n in range(M_array[l, k]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :k]) + n
                    filestring += f"### OUTER LEVEL {l-1}, INNER LEVEL {k-1}, SAMPLE {sample_number}\n\n"
                    filestring += f"cp {fname}_olevel{l-1}_ilevel{k-1}_sample{sample_number}.swn INPUT\n"
                    filestring += "$swan_omp_exe\n"
                    filestring += f"cp PRINT {fname}_olevel{l-1}_ilevel{k-1}_sample{sample_number}.prt\n"
                    filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
                for n in range(M_array[l, k]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :k]) + n
                    filestring += f"### OUTER LEVEL {l}, INNER LEVEL {k-1}, SAMPLE {sample_number}\n\n"
                    filestring += f"cp {fname}_olevel{l}_ilevel{k-1}_sample{sample_number}.swn INPUT\n"
                    filestring += "$swan_omp_exe\n"
                    filestring += f"cp PRINT {fname}_olevel{l}_ilevel{k-1}_sample{sample_number}.prt\n"
                    filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
                for n in range(M_array[l, k]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :k]) + n
                    filestring += f"### OUTER LEVEL {l-1}, INNER LEVEL {k}, SAMPLE {sample_number}\n\n"
                    filestring += f"cp {fname}_olevel{l-1}_ilevel{k}_sample{sample_number}.swn INPUT\n"
                    filestring += "$swan_omp_exe\n"
                    filestring += f"cp PRINT {fname}_olevel{l-1}_ilevel{k}_sample{sample_number}.prt\n"
                    filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
                for n in range(M_array[l, k]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :k]) + n
                    filestring += f"### OUTER LEVEL {l}, INNER LEVEL {k}, SAMPLE {sample_number}\n\n"
                    filestring += f"cp {fname}_olevel{l}_ilevel{k}_sample{sample_number}.swn INPUT\n"
                    filestring += "$swan_omp_exe\n"
                    filestring += f"cp PRINT {fname}_olevel{l}_ilevel{k}_sample{sample_number}.prt\n"
                    filestring += "rm -f INPUT\nrm -f PRINT\nrm -f norm_end\n\n"
        
        f = open(path + fname + '.sh', 'x')
        f.write(filestring)


## DOWNLOAD AND MANIPULATE ECMWF DATA ##

class ecmwfdata:
    """Contains tools to download, write, and manipulate data from ECMWF. Also contains attributes of compatible wind input for SWAN.

    - correct_wind_time_attrs (dict):   attributes (metadata) of time dimension for SWAN wind input;
    - correct_wind_x_attrs (dict):      attributes (metadata) of longitude dimension for SWAN wind input;
    - correct_wind_y_attrs (dict):      attributes (metadata) of latitude dimension for SWAN wind input;
    - correct_wind_z_attrs (dict):      attributes (metadata) of 'mean height above sea level' data variable for SWAN wind input;
    - correct_wind_crs_attrs (dict):    attributes (metadata) of coordinate reference system data variable for SWAN input;
    - correct_wind_u_attrs (dict):      attributes (metadata) of eastward wind for SWAN wind input;
    - correct_wind_v_attrs (dict):      attributes (metadata) of northward wind for SWAN wind input;
    - correct_wind_attrs (dict):        attributes (metadata) of SWAN wind input dataset as a whole;
    
    - generate_ecmwfapirc:              generates a .ecmwfapirc-file, to "log-in" to the ECMWF data servers with your account;  
    - retrieve_tigge_wind_data:         downloads TIGGE-data of a selection of dates, ensemble members and geographical area from the ECMWF data server;
    - tigge_request:                    downloads a TIGGE-ensemble forecast of a specific date;
    - extract_wind_ensemble_from_tigge: converts a TIGGE-GRIB-file into a netcdf-file for every ensemble member; makes it compatible with SWAN;
    - convert_era5wind_to_swaninput:    converts ERA5 wind data to a file compatible with SWAN;
    - tpar_from_era5wavefield:          writes *.tpar-boundary data files from an ERA5 wave dataset;
    
    """

    # attributes of a wind input file that is compatible with SWAN

    correct_wind_time_attrs = {'long_name': 'time', 'standard_name': 'time', 'axis': 'T'}
    correct_wind_x_attrs = {'units': 'degrees_east', 'long_name': 'x coordinate according to WGS 1984', 'standard_name': 'longitude', 'axis': 'X'}
    correct_wind_y_attrs = {'units': 'degrees_north', 'long_name': 'y coordinate according to WGS 1984', 'standard_name': 'latitude', 'axis': 'Y'}
    correct_wind_z_attrs = {'units': 'meters', 'long_name': 'height above mean sea level'}
    correct_wind_crs_attrs = {
        'units': 'degrees_north',
        'long_name': 'coordinate reference system',
        'grid_mapping_name': 'latitude_longitude',
        'longitude_of_prime_meridian': 0.0,
        'semi_major_axis': 6378137.0,
        'inverse_flattening': 298.257223563,
        'epsg_code': 'EPSG:4326',
        'proj4_params': '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
        'crs_wkt': 'GEOGCS["WGS 84",\n        DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],\n        PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],\n        UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],\n        AUTHORITY["EPSG","4326"]] '
    }
    correct_wind_u_attrs = {'units': 'm/s', 'long_name': 'eastward_wind', 'standard_name': 'eastward_wind', 'grid_mapping': 'crs'}
    correct_wind_v_attrs = {'units': 'm/s', 'long_name': 'northward_wind', 'standard_name': 'northward_wind', 'grid_mapping': 'crs'}
    correct_wind_attrs = {'Conventions': 'CF-1.6,UGRID-0.9', 'Metadata_Conventions': 'Unidata Dataset Discovery v1.0'}


    @staticmethod
    def generate_ecmwfapirc_file(path, key, email):
        """Generate a .ecmwfapirc-file to log-in to the ECMWF data server. Once you have made a free account on the ECMWF website,
        you can find you API key there. This contains the information this function requires as input. Recommended to save this file to
        C:\\Users\\{username}\\.

        - path (str):   path to which the .ecmwfapirc-file will be written (recommended 'C:\\Users\\{username}\\');
        - key (str):    personal ECMWF API key;
        - email (str):  registered email at ECMWF.
        """
        url_api = "https://api.ecmwf.int/v1"

        # construct string
        filestring = '{\n    "url"   : "' + url_api + '",\n    "key"   : "' + key + '",\n    "email" : "' + email + '"\n}' 
        print(filestring)
        f = open(path + '.ecmwfapirc', 'x')
        f.write(filestring)


    @staticmethod
    def retrieve_tigge_wind_data(dates, members, hours, min_lon, max_lon, min_lat, max_lat):
        """Downloads TIGGE-data of perturbed forecasts from a set of dates.
        
        - dates (list):     list of dates from which to download forecasts (yyyy-mm-dd);
        - members (list):   list of indices indicating which ensemble members should be included;
        - hours (list):     list indicating how long the forecast should be per forecast data;
        - min_lon (float):  west boundary of domain (multiple of 0.5 degrees);
        - max_lon (float):  east boundary of domain (multiple of 0.5 degrees);
        - min_lat (float):  south boundary of domain (multiple of 0.5 degrees);
        - max_lat (float):  north boundary of domain (multiple of 0.5 degrees).        
        """
        server = ECMWFDataServer()

        # datestring = ""
        # for date in dates:
        #     datestring += date + '/'
        # datestring = datestring[:-1]
        
        for i, date in enumerate(dates):
            target = "ecmwf_sfc_%s_00.grb" % (date)
            ecmwfdata.tigge_request(date, members, hours[i], server, target, min_lon, max_lon, min_lat, max_lat)


    @staticmethod
    def tigge_request(date, members, hour, server, target, min_lon, max_lon, min_lat, max_lat):
        """Downloads TIGGE-data of perturbed forecasts from a specific date from the ECMWF data server.
        
        - date (str):       date from which to download forecast (yyyy-mm-dd);
        - members (list):   list of indices indicating which ensemble members should be included;
        - hour (int):       duration of the forecast; should be a multiple of 6;
        - server:           ECMWFDataServer object initialised in e.g. retrieve_tigge_wind_data;
        - target (str):     name of the grib-file that should be downloaded;
        - min_lon (float):  west boundary of domain (multiple of 0.5 degrees);
        - max_lon (float):  east boundary of domain (multiple of 0.5 degrees);
        - min_lat (float):  south boundary of domain (multiple of 0.5 degrees);
        - max_lat (float):  north boundary of domain (multiple of 0.5 degrees);
        """

        memberstring = ''
        for member in members:
            memberstring += str(member) + '/'
        memberstring = memberstring[:-1]

        area = f'{max_lat}/{min_lon}/{min_lat}/{max_lon}' # north - west - south - east boundary

        hourstring = ''
        for i in range(hour // 6 + 1):
            hourstring += str(6*i) + '/'
        hourstring = hourstring[:-1] # delete final '/'
        

        server.retrieve({
            "class": "ti",
            "dataset": "tigge",
            "date": date,
            "expver": "prod",
            "grid": "0.5/0.5",
            "levtype": "sfc",
            "number": memberstring,
            "origin": "ecmwf",
            "param": "165/166", #u- and v-components of 10m wind
            "area": area,
            "step": hourstring, 
            "target": target,
            "time": "00:00:00",
            "type": "pf",
        }) 


    @staticmethod
    def extract_wind_ensemble_from_tigge(grib_path, grib_name, ensemble_member_path, ensemble_member_name):
        """Given a GRIB-file of TIGGE-data, opens and writes this data to a netcdf-file per ensemble member. The function will make these netcdf-files usable
        as input files for the SWAN model.
        
        - grib_path (str):              path to the GRIB-file;
        - grib_name (str):              name of the GRIB-file (with extension);
        - ensemble_member_path (str):   path to the folder where ensemble member data should be written to;
        - ensemble_member_name (str):   name of the ensemble member files (without ensemble member index and extension; function adds these automatically).
        """
            
        # Open GRIB-data
        print('Opening Data...\n\n')
        full_data = xr.open_dataset(grib_path + grib_name, engine='cfgrib')

        start_time = full_data['time'] # this is the start time as an np.datetime64-object
        time_steps = full_data['step'] # this are the time steps as np.timedelta64-objects

        # Metadata for SWAN-input encoding
        time_unit  = "minutes since 1970-01-01 00:00:00.0 +0000"
        fill_value = 9.96921e+36

        encoding_SWAN = {'time': {'dtype': 'float64', '_FillValue': fill_value, 'units': time_unit},
                     'x':    {'dtype': 'float64', '_FillValue': fill_value},
                     'y':    {'dtype': 'float64', '_FillValue': fill_value},
                     'z':    {'dtype': 'float64', '_FillValue': fill_value},
                     'crs':  {'dtype': 'int32'},
                     'eastward_wind':  {'dtype': 'float32', '_FillValue': -999.},
                     'northward_wind': {'dtype': 'float32', '_FillValue': -999.},
                     }

        # Construct coordinates

        print('Constructing coordinate arrays...\n\n')
        new_time = xr.DataArray(full_data['valid_time'].values, dims = ('time'), name='time')
        x = xr.DataArray(full_data['longitude'].values, dims=('x'), name='x')
        y = xr.DataArray(full_data['latitude'].values, dims=('y'), name='y')

        new_time.attrs = ecmwfdata.correct_wind_time_attrs
        x.attrs = ecmwfdata.correct_wind_x_attrs
        y.attrs = ecmwfdata.correct_wind_y_attrs

        # Data variables
        print('Constructing data variables...\n\n')
        z_field = np.empty((full_data['latitude'].values.shape[0], full_data['longitude'].values.shape[0]))
        z_field[:] = 10.
        z = xr.DataArray(z_field, dims=('y', 'x'), coords={'y': y, 'x': x}, name='z')
        z.attrs = ecmwfdata.correct_wind_z_attrs

        crs = xr.DataArray(0, name='crs')
        crs.attrs = ecmwfdata.correct_wind_crs_attrs
        # Wind data
        for n in range(full_data['number'].values.shape[0]):
            eastward_wind = xr.DataArray(full_data['u10'].values[n, :, :, :], dims=('time', 'y', 'x'), coords={'time': new_time, 'y': y, 'x': x}, name='eastward_wind')
            eastward_wind.attrs = ecmwfdata.correct_wind_u_attrs

            northward_wind = xr.DataArray(full_data['v10'].values[n, :, :, :], dims=('time', 'y', 'x'), coords={'time': new_time, 'y': y, 'x': x}, name='eastward_wind')
            northward_wind.attrs = ecmwfdata.correct_wind_v_attrs

            wind_field = xr.Dataset({'z': z, 'crs': crs, 'eastward_wind': eastward_wind, 'northward_wind': northward_wind})
            wind_field.attrs = ecmwfdata.correct_wind_attrs
            wind_field = wind_field[['time', 'y', 'x', 'z', 'crs', 'eastward_wind', 'northward_wind']]

            wind_field.to_netcdf(ensemble_member_path + ensemble_member_name + str(n) + '.nc', format='NETCDF3_CLASSIC', encoding=encoding_SWAN)


    @staticmethod
    def convert_era5wind_to_swaninput(path):
        """Converts ERA5 wind files (*.nc) to wind input files (*.nc) that are compatible with SWAN. Strongly based on code by Stendert Laan 
        (stendert.laan@deltares.nl). This function works with ERA5 wind files that do not have dimension 'number', so ensemble members should
        already have been extracted. The new wind input will be written to the same path as the ERA5_data, with name 
        '{name_of_era5}_compatible.nc'.
        
        - path (str):       path of the ERA5 dataset.
        
        """

        correct_wind_data = xr.open_dataset('P:\\1230882-emodnet_hrsm\\vanCas\\vancg\\boundary_conditions\\meteo\\anonymus_hwang_wind.nc')
        era5_wind_field = xr.open_dataset(path, engine='netcdf4')

        time_unit  = "minutes since 1970-01-01 00:00:00.0 +0000"
        fill_value = 9.96921e+36

        # Coordinates

        time = xr.DataArray(era5_wind_field['time'].values, dims=('time'), name='time')
        x = xr.DataArray(era5_wind_field['longitude'].values, dims=('x'), name='x')
        y = xr.DataArray(era5_wind_field['latitude'].values, dims=('y'), name='y')

        time.attrs = correct_wind_data['time'].attrs
        x.attrs = correct_wind_data['x'].attrs
        y.attrs = correct_wind_data['y'].attrs

        # Data variables

        z_field = np.empty((era5_wind_field['latitude'].values.shape[0], era5_wind_field['longitude'].values.shape[0]))
        z_field[:] = 10.
        z = xr.DataArray(z_field, dims=('y', 'x'), coords={'y': y, 'x': x}, name='z')
        z.attrs = correct_wind_data['z'].attrs

        crs = xr.DataArray(0, name='crs')
        crs.attrs = correct_wind_data['crs'].attrs

        eastward_wind = xr.DataArray(era5_wind_field['u10'].values, dims=('time', 'y', 'x'), coords={'time': time, 'y': y, 'x': x}, name='eastward_wind')
        eastward_wind.attrs = correct_wind_data['eastward_wind'].attrs

        northward_wind = xr.DataArray(era5_wind_field['v10'].values, dims=('time', 'y', 'x'), coords={'time': time, 'y': y, 'x': x}, name='northward_wind')
        northward_wind.attrs = correct_wind_data['northward_wind'].attrs

        swn_wind_field = xr.Dataset({'z': z, 'crs': crs, 'eastward_wind': eastward_wind, 'northward_wind': northward_wind})
        swn_wind_field.attrs['Conventions']             = 'CF-1.6,UGRID-0.9'
        swn_wind_field.attrs['Metadata_Conventions']    = 'Unidata Dataset Discovery v1.0'
        swn_wind_field = swn_wind_field[['time', 'y', 'x', 'z', 'crs', 'eastward_wind', 'northward_wind']]
        
        encoding_SWAN = {'time': {'dtype': 'float64', '_FillValue': fill_value, 'units': time_unit},
                        'x':    {'dtype': 'float64', '_FillValue': fill_value},
                        'y':    {'dtype': 'float64', '_FillValue': fill_value},
                        'z':    {'dtype': 'float64', '_FillValue': fill_value},
                        'crs':  {'dtype': 'int32'},
                        'eastward_wind':  {'dtype': 'float32', '_FillValue': -999.},
                        'northward_wind': {'dtype': 'float32', '_FillValue': -999.},
                        }

        swn_wind_field.to_netcdf(path[:-3]+'_compatible.nc', format='NETCDF3_CLASSIC', encoding=encoding_SWAN)

    @staticmethod
    def tpar_from_era5wavefield(era5wavefield, tparpath, tparname):
        """Converts an ERA5 reanalysis dataset into nonstationary boundary condition files (*.tpar) for use as SWAN input. Geographical
        domain of SWAN computation and ERA5 dataset must be the same. The ERA5 dataset must contain (at least) the parameters 'Significant height
        of combined wind waves and swell', 'Significant height of first swell partition', 'Significant height of second swell partition',
        'Significant height of third swell partition', 'Significant height of wind waves', 'Peak wave period', 'Mean direction of wind waves', 
        'Mean wave direction of first swell partition', 'Mean wave direction of second swell partition', 'Mean wave direction of third swell
        partition' and 'Wave Spectral Directional width'.
        
        - era5wavefield (xArray dataset):   dataset containing all wave parameters;
        - tparpath (str):                   path to which the boundary condition files will be written;
        - tparname (str):                   file name of the *.tpar files (without location index or extension).
        """

        lon = era5wavefield['longitude']
        lat = era5wavefield['latitude']
        time = era5wavefield['time']
        
        hstotal = era5wavefield['swh']
        hswind = era5wavefield['shww']
        hswe1 = era5wavefield['p140121']
        hswe2 = era5wavefield['p140124']
        hswe3 = era5wavefield['p140127']

        per = era5wavefield['pp1d']

        dirwind = era5wavefield['mdww']
        dirswe1 = era5wavefield['p140122']
        dirswe2 = era5wavefield['p140125']
        dirswe3 = era5wavefield['p140128']

        dspr = era5wavefield['wdw'] * 180/(2*np.sqrt(np.pi)) # Conversion factor to something SWAN can interpret
        # Suggestion from Caroline Gautier: dspr.values = np.maximum(np.full(dspr.values.shape, 31.5), dspr.values)

        # Determine peak wave direction by comparing significant wave heights of wind waves and all swell partitions
        maximum_hs = np.array([hswind.values, hswe1.values, hswe2.values, hswe3.values]).argmax(axis=0) # indicator which wave type is maximal

        peak_dir_values = np.zeros(hswind.values.shape)
        peak_dir_values = np.where(maximum_hs != 0, peak_dir_values, dirwind)
        peak_dir_values = np.where(maximum_hs != 1, peak_dir_values, dirswe1) 
        peak_dir_values = np.where(maximum_hs != 2, peak_dir_values, dirswe2)
        peak_dir_values = np.where(maximum_hs != 3, peak_dir_values, dirswe3)

        ## WRITE *.TPAR-FILES ##

        # West boundary
        for i in range(lat.values.shape[0]):
            if np.isnan(hstotal[0, i, 0]) or np.isnan(per[0,i,0]) or np.isnan(dspr[0,i,0]) or np.isnan(peak_dir_values[0,i,0]): # Don't create a *.tpar-file if current pixel is a dry pixel
                continue

            string = "TPAR\n"
            for n, t in enumerate(time.values):
                iso = ppt.string.datetime_to_isostring(t)

                hs_point = float(hstotal.values[n, i, 0])
                per_point = float(per.values[n, i, 0])
                pdir_point = float(peak_dir_values[n, i, 0])
                dspr_point = float(dspr.values[n, i, 0])
                string += iso + f" {np.round(hs_point, 1)} {np.round(per_point, 1)} {np.round(pdir_point, 0)} {np.round(dspr_point,0)}\n"
                # rounding is necessary because SWAN doesn't work if too many decimals are present
            try:
                f = open(tparpath + f'{tparname}_w_{i}.tpar', 'x')
            except FileExistsError:
                print('Warning: overwriting a file')
                f = open(tparpath + f'{tparname}_w_{i}.tpar', 'w')
            f.write(string)

        # North boundary
        
        for i in range(lon.values.shape[0]):
            if np.isnan(hstotal[0, 0, i])or np.isnan(per[0,0,i]) or np.isnan(dspr[0,0,i]) or np.isnan(peak_dir_values[0,0,i]): # Don't create a *.tpar-file if current pixel is a dry pixel
                continue

            string = "TPAR\n"
            for n, t in enumerate(time.values):
                iso = ppt.string.datetime_to_isostring(t)

                hs_point = float(hstotal.values[n, 0, i])
                per_point = float(per.values[n, 0, i])
                pdir_point = float(peak_dir_values[n, 0, i])
                dspr_point = float(dspr.values[n, 0, i])
                string += iso + f" {np.round(hs_point, 1)} {np.round(per_point, 1)} {np.round(pdir_point, 0)} {np.round(dspr_point,0)}\n"
            try:
                f = open(tparpath + f'{tparname}_n_{i}.tpar', 'x')
            except FileExistsError:
                print('Warning: overwriting a file')
                f = open(tparpath + f'{tparname}_n_{i}.tpar', 'w')
            f.write(string)

        # East boundary
       
        for i in range(lat.values.shape[0]):
            if np.isnan(hstotal[0, i, -1]) or np.isnan(per[0,i,-1]) or np.isnan(dspr[0,i,-1]) or np.isnan(peak_dir_values[0,i,-1]): # Don't create a *.tpar-file if current pixel is a dry pixel
                continue

            string = "TPAR\n"
            for n, t in enumerate(time.values):
                iso = ppt.string.datetime_to_isostring(t)

                hs_point = float(hstotal.values[n, i, -1])
                per_point = float(per.values[n, i, -1])
                pdir_point = float(peak_dir_values[n, i, -1])
                dspr_point = float(dspr.values[n, i, -1])
                string += iso + f" {np.round(hs_point, 1)} {np.round(per_point, 1)} {np.round(pdir_point, 0)} {np.round(dspr_point,0)}\n"
            try:
                f = open(tparpath + f'{tparname}_e_{i}.tpar', 'x')
            except FileExistsError:
                print('Warning: overwriting a file')
                f = open(tparpath + f'{tparname}_e_{i}.tpar', 'w')
            f.write(string)

        # South boundary
        
        for i in range(lon.values.shape[0]):
            if np.isnan(hstotal[0, -1, i]) or np.isnan(per[0,-1,i]) or np.isnan(dspr[0,-1,i]) or np.isnan(peak_dir_values[0,-1,i]): # Don't create a *.tpar-file if current pixel is a dry pixel
                continue

            string = "TPAR\n"
            for n, t in enumerate(time.values):
                iso = ppt.string.datetime_to_isostring(t)

                hs_point = float(hstotal.values[n, -1, i])
                per_point = float(per.values[n, -1, i])
                pdir_point = float(peak_dir_values[n, -1, i])
                dspr_point = float(dspr.values[n, -1, i])
                string += iso + f" {np.round(hs_point, 1)} {np.round(per_point, 1)} {np.round(pdir_point, 0)} {np.round(dspr_point,0)}\n"
            try:
                f = open(tparpath + f'{tparname}_s_{i}.tpar', 'x')
            except FileExistsError:
                print('Warning: overwriting a file')
                f = open(tparpath + f'{tparname}_s_{i}.tpar', 'w')
            f.write(string)

    
##############
## OLD CODE ##
##############

# def extract_ensemble_members(file, fname, path):
#     """Extracts the ensemble members from an ERA5-ensemble member dataset (in netcdf-format). Assumes coordinates to be:
#     'time', 'number', 'longitude', 'latitude'.
    
#     - file (str):       file path of the ERA5-ensemble member dataset;
#     - fname (str):      naming convention of the individual ensemble members; eventual file names will be '{fname}_{number}.nc';
#     - path (str):       path to which the individual ensemble members should be written.
    
#     """
#     f = xr.open_dataset(file, engine='netcdf4')
#     ensemble_size = f['number'].values.shape[0]
#     for n in range(ensemble_size):
#         data_vars_dict = dict()
#         for var_name in f.keys():
#             data_vars_dict[var_name] = (['time', 'latitude', 'longitude'], f[var_name][:, n, :, :].values)
#         coords = {'time': ('time', f['time'].values),
#                   'longitude': ('longitude', f['longitude'].values),
#                   'latitude': ('latitude', f['latitude'].values)}
#         attrs = {'Conventions': 'CF-1.6'}
#         data = xr.Dataset(data_vars=data_vars_dict, coords=coords, attrs=attrs)
#         data.to_netcdf(path+fname+f'_{n}.nc'
    
# @staticmethod
# def write_MLMC_EPS_input_files(input_path, fname, wind_path_and_wind_filename, num_samples, N_list, approx='spatial_resolution3',
#                            tpar_path=None):
#     # Check if N_list is valid

#     if sum(N_list) > num_samples:
#         raise ValueError("Not enough samples available. Decrease numbers in N_list.")
#     elif sum(N_list) < num_samples:
#         not_all_samples_continue = input("Warning: not all available samples are used. Continue anyway? [y/n] ")
#         if not_all_samples_continue == 'y' or not_all_samples_continue == 'Y' or not_all_samples_continue == 'yes' or not_all_samples_continue == 'Yes':
#             num_samples = sum(N_list)
#         else:
#             return
        
#     ## Approximation using full resolution, half resolution and quarter resolution ##################################

#     if approx == 'spatial_resolution3':

#         # Check if N_list is valid for 'spatial_resolution3'-approximation type

#         if len(N_list) != 3:
#             raise ValueError("Length of N_list should be 3 when using the 'spatial_resolution3'-approximations.")

#         # Input options
#         options = input('Default parameters or input manually? [default/manual] ')
#         while options != 'default':
#             print('Sorry, I did not hear you correctly. Did you mean default? ')
#             options = input('Default parameters or input manually? [default/manual] ')
            
#         if options == 'default':
#             start_iso = '20131205.000000'
#             end_iso = '20131207.210000'
#             timestep = 1
#             project = 'SWAN-EPS'
#             freq_lower_bound = .03
#             boundary_direction_type = 'PEAK'
#             generation = 3
#             npnts = 98
#             maxitns = 20
#         else:
#             # not coded yet
#             return
        
#         run_counter = 0

#         # Lowest approximation level

#         for n in range(N_list[0]):
#             sample_number = n
#             swnstr = Swanstring(f'{fname}_level0_{sample_number}', project, ppt.int_to_threedigitstring(run_counter), start_iso,
#                                 end_iso, timestep)
#             swnstr.input(105, 120, '../../geometry/swan-ns-j22_6-v1a_adjust.bot', windfile=f'../../boundary_conditions/meteo/members/sinterklaasstorm_wind_ensemble_{sample_number}.nc',
#                         spec_lowerbound=freq_lower_bound)
#             # Boundary conditions; can only write this code if I know how it works

#             swnstr.physics(generation=generation)
#             swnstr.numerical_parameters(npnts=npnts, maxitns=maxitns)
#             swnstr.output(f'../../geometry/output_locations/', f'../../output/0/', f'eps0_level0_sample{sample_number}')
#             swnstr.write_to_file("P:\\1230882-emodnet_hrsm\\vanCas\\Ensemble_Forecast\\input\\0")

#         # Higher approximation levels

#         for l in range(1, 3):
#             for n in range(N_list[l]):
#                 sample_number = sum(N_list[:l]) + n

#                 swnstr0 = Swanstring(f'{fname}_level{l-1}_{sample_number}', project, ppt.int_to_threedigitstring(run_counter), start_iso,
#                                 end_iso, timestep)
#                 swnstr0.input(105, 120, '../../geometry/swan-ns-j22_6-v1a_adjust.bot', windfile=f'../../boundary_conditions/meteo/members/sinterklaasstorm_wind_ensemble_{sample_number}.nc',
#                             spec_lowerbound=freq_lower_bound)
#                 # Boundary conditions; can only write this code if I know how it works

#                 swnstr0.physics(generation=generation)
#                 swnstr0.numerical_parameters(npnts=npnts, maxitns=maxitns)
#                 swnstr0.output(f'../../geometry/output_locations/', f'../../output/0/', f'eps0_level{l-1}_sample{sample_number}')
#                 swnstr0.write_to_file("P:\\1230882-emodnet_hrsm\\vanCas\\Ensemble_Forecast\\input\\0")

#                 swnstr1 = Swanstring(f'{fname}_level{l}_{sample_number}', project, ppt.int_to_threedigitstring(run_counter), start_iso,
#                                 end_iso, timestep)
#                 swnstr1.input(105, 120, '../../geometry/swan-ns-j22_6-v1a_adjust.bot', windfile=f'../../boundary_conditions/meteo/members/sinterklaasstorm_wind_ensemble_{sample_number}.nc',
#                             spec_lowerbound=freq_lower_bound)
#                 # Boundary conditions; can only write this code if I know how it works

#                 swnstr1.physics(generation=generation)
#                 swnstr1.numerical_parameters(npnts=npnts, maxitns=maxitns)
#                 swnstr1.output(f'../../geometry/output_locations/', f'../../output/0/', f'eps0_level{l}_sample{sample_number}')
#                 swnstr1.write_to_file("P:\\1230882-emodnet_hrsm\\vanCas\\Ensemble_Forecast\\input\\0")

