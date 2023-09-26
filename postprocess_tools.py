import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import copy
from scipy.interpolate import RegularGridInterpolator
from skimage.restoration import inpaint


## STANDARD STATISTICAL FUNCTIONS ##


def RMSE(arr1, arr2):
    effective_size = np.count_nonzero(~np.isnan((arr1-arr2).values))
    square_sum = np.nansum(np.power((arr1.values-arr2.values),2))
    return np.sqrt(square_sum/effective_size)

## PLOTTING TOOLS ##

class plot:
    """Contains methods to visualise SWAN output data. For most methods in this class, a figure and axis object have to be predefined. In this case,
    the methods of the class add plots to one of the axis objects of the figure. Only the method 'plot_spectrum2d' does not rely on a predefined axis object.
    It also follows that 'plt.show()' has to be added for most of these visualisations to show up.
    
    - plot_colormap:            plot a colorvalued map of a quantity varying in latitude and longitude;
    - plot_point_on_colormap:   given an axis object with a colormap (such as the output of 'plot_colormap'), plots a point with a given latitude and longitude on the map;
    - plot_arrows_on_colormap:  given an axis object with a colormap, plots quivers on that map (either vectors or angles);
    - plot_timeseries:          plot a time series of a quantity at a certain geographical location;
    - plot_spectrum1d:          plot the 1d-energy density spectrum at a geographical location and at some given timestep;
    - plot_spectrum2d:          plot the 2d-energy density spectrum at a geographical location and at some given timestep.
    """

    @staticmethod
    def plot_colormap(fig, ax, longitude, latitude, valmap, title, quantitystr=None, vmin=None, vmax=None, cmap='viridis'):
        """Plots a colormap on a given axis object in a given matplolib.pyplot figure.

        - fig (plt.Figure):                             figure in which the colormap will be plotted;
        - ax (plt.Axes):                                axis on which the colormap will be plotted;
        - longitude (xArray DataArray or np.ndarray):   longitude coordinates;
        - latitude (xArray DataArray or np.ndarray):    latitude coordinates;
        - valmap (xArray DataArray or np.ndarray):      values that should be plotted in the colormap;
        - title (str):                                  title of the plot;
        - quantitystr (str):                            name of the quantity being plotted (with unit);
        - vmin (float):                                 lowest value that is assigned a color;
        - vmax (float):                                 highest value that is assigned a color;
        - cmap (str):                                   name of the colormap that will be used.
        """
        im = ax.pcolormesh(longitude, latitude, valmap, cmap=cmap, zorder=1, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(im, pad=.1)

        if quantitystr:
            cbar.set_label(quantitystr)

        ax.set_title(title)


    @staticmethod
    def plot_point_on_colormap(ax, poi_lon, poi_lat, color):
        """Plots a point on an axis object with a colormap.
        
        - ax (plt.Axes):    axis containing the colormap;
        - poi_lon (float):  longitude coordinate of point;
        - poi_lat (float):  latitude coordinate of point;
        - color (str):      color of the point.
        """
        poi = ax.scatter(poi_lon, poi_lat, color=color, marker='o', zorder=4)


    @staticmethod
    def plot_arrows_on_colormap(ax, lon, lat, stride, quantity, color='black', mode='angle'):
        """Plots arrows representing e.g. wave direction or wind direction on a colormap. Input is either an array of angles or vectors.
        
        - ax (plt.Axes):                                    axis containing the colormap;
        - lon (xArray DataArray or np.ndarray):             longitude coordinates;
        - lat (xArray DataArray or np.ndarray):             latitude coordinates;
        - stride (int):                                     number of grid cells between each arrow;
        - quantity (xArray DataArray, np.ndarray or tuple): data that defines the direction of the arrow. If mode == 'angle', this should be an array (xArray or numpy) and
                                                            if mode == 'vector', this should be a tuple of arrays (X, Y) containing the x- and y-coordinates respectively;
        - color (str):                                      color of the arrows;
        - mode (str):                                       'angle' or 'vector'.
        """
        lon_quiver, lat_quiver = np.meshgrid(lon[::stride], lat[::stride])
        if mode == 'angle':
            angle_rad = np.radians(quantity[::stride, ::stride])

            quiv = ax.quiver(lon_quiver, lat_quiver, np.cos(angle_rad), np.sin(angle_rad), pivot='middle', color=color, 
                            zorder=3)
        elif mode == 'vector':
            magnitude = np.sqrt(np.power(quantity[0].values, 2) + np.power(quantity[1].values, 2))
            norm_x = quantity[0].values / magnitude
            norm_y = quantity[1].values / magnitude

            quiv = ax.quiver(lon_quiver, lat_quiver, norm_x, norm_y, pivot='middle', color=color, zorder = 3)
        

    @staticmethod
    def plot_timeseries(ax, lon, lat, poi_lon, poi_lat, qoi, time, title=None, quantitystr=None, label=None, color=None):
        """Plots a time series on a given axis object.
        
        - ax (plt.Axes):                            axis containing the time series;
        - lon (xArray DataArray):                   longitude coordinates;
        - lat (xArray DataArray):                   latitude coordinates;
        - poi_lon (float):                          longitude coordinate of the location of interest;
        - poi_lat (float):                          latitude coordinate of the location of interest;
        - qoi (xArray DataArray or np.ndarray):     quantity that will be plotted; 
        - time (xArray DataArray or np.ndarray):    time coordinates;     
        - title (str):                              title of the time series;
        - quantitystr (str):                        name of the quantity of interest (with unit);
        - label (str):                              label of this time series; relevant if multiple time series are plotted;
        - color (str):                              color of the line;          
        """
        lon_index = np.argmin(np.absolute(lon.values - np.full(lon.shape, poi_lon)))
        lat_index = np.argmin(np.absolute(lat.values - np.full(lat.shape, poi_lat)))
        if label:
            ax.plot(time, qoi[:,lon_index, lat_index], label=label, color=color)
        else:
            ax.plot(time, qoi[:,lon_index, lat_index], color=color)
        ax.set_xlabel('Time')
        if quantitystr:
            ax.set_ylabel(quantitystr)
        if title:
            ax.set_title(title)

        
    @staticmethod
    def plot_spectrum1d(ax, frequency, spectrum, toi_index, poi_index, label=None, title=None, color='blue'):
        """Plots the one-dimensional energy density spectrum at a certain point in time and space on a given axis object.
        
        - ax (plt.Axes):                                axis containing the spectrum plot;
        - frequency (xArray DataArray or np.ndarray):   frequency coordinates;
        - spectrum (xArray DataArray or np.ndarray):    energy density depending on time, space and frequency;
        - toi_index (int):                              time index of interest;
        - poi_index (int):                              index of the point of interest in the list of output locations;
        - label (str):                                  label of this energy density spectrum plot;
        - title (str):                                  title of this energy density spectrum plot;
        - color (str):                                  color of the line.
        """
        if label:
            ax.plot(frequency, spectrum[toi_index, poi_index, :], label=label, color=color)
        else:
            ax.plot(frequency, spectrum[toi_index, poi_index, :], color=color)
        
        if title:
            ax.set_title(title)
        
        ax.set_xlabel(f'Frequency [Hz]')
        ax.set_ylabel(f'Energy density [m^2 s]')


    @staticmethod
    def plot_spectrum2d(frequency, direction, spectrum, toi_index, poi_index, title=None, savestr=None, cmap='viridis'):
        """Generates a polar plot of the two-dimensional energy density spectrum at a certain point in time and space. This method
        takes care of defining the figure and axis objects and 'plt.show()'. Directional convention is nautical.
        
        - frequency (xArray DataArray or np.ndarray):   frequency coordinates;
        - direction (xArray DataArray):                 direction coordinates;
        - spectrum (xArray DataArray):                  energy density values;
        - toi_index (int):                              time index of interest;
        - poi_index (int):                              index of the point of interest in the list of output locations;
        - title (str):                                  title of this plot;
        - savestr (str):                                path and filename to which the figure will be saved. If None, the figure will not be saved;
        - cmap (str):                                   name of the colormap.
        """
        spectrum_data = spectrum.isel(time=toi_index, points=poi_index).values
        # Convert azimuthal direction to northing direction (clockwise from north, where waves propagate towards in stead of from)
        direction_degrees = np.degrees(direction.values)
        direction_northing = 360 - direction_degrees
        direction_northing_rad = np.radians(direction_northing)

        figs, axs = plt.subplots(subplot_kw={'projection': 'polar'})
        im = axs.pcolormesh(direction_northing_rad, frequency, spectrum_data, cmap=cmap)
        axs.set_theta_zero_location('N')

        cbar = figs.colorbar(im, pad=.1)
        cbar.set_label('Directional energy density [m^2 s]')
        if title:
            axs.set_title(title)
        axs.set_xlabel('Wave direction [degrees]')
        axs.set_ylabel('Frequency [Hz]')
        axs.yaxis.labelpad=27
        axs.set_rlabel_position(90)
        
        if savestr:
            figs.savefig(savestr)

        plt.show()


## STRING MANIPULATION TOOLS ##

class string:
    """Convenient methods to manipulate strings specific to SWAN-application.
    
    - datetime_to_isostring:        converts np.datetime64 object into a string with the ISO-notation;
    - isostring_to_datestring:      converts ISO-notation into a readable date text;
    - num_to_month:                 converts month number into month name;
    - int_to_three_digit_string:    converts integer (0-999) into a string consisting of three characters.
    
    """
    @staticmethod
    def datetime_to_isostring(time):
        string = str(time)
        return string[:4] + string[5:7] + string[8:10] + '.' + string[11:13] + string[14:16] + string[17:19]

    @staticmethod
    def isostring_to_datestring(iso, year=False):
        day = iso[6:8]
        time = iso[9:11] + '.' + iso[11:13]
        month = iso[4:6]
        monthstr = string.num_to_month(month)

        if year:
            yearstr = iso[:4]
            return f'{monthstr} {day} {yearstr}, {time}'
        else:
            return f'{monthstr} {day}, {time}'
        
    @staticmethod
    def num_to_month(num):
        if num == '01':
            return 'January'
        elif num == '02':
            return 'February'
        elif num == '03':
            return 'March'
        elif num == '04':
            return 'April'
        elif num == '05':
            return 'May'
        elif num == '06':
            return 'June'
        elif num == '07':
            return 'July'
        elif num == '08':
            return 'August'
        elif num == '09':
            return 'September'
        elif num == '10':
            return 'October'
        elif num == '11':
            return 'November'
        elif num == '12':
            return 'December'
        else:
            raise ValueError(f'Could not convert {num} to a month.')
    

    @staticmethod
    def int_to_threedigitstring(num):
        return '0' * (3-len(str(num))) + str(num)


## PREDICT MONTE CARLO METHOD RUNTIMES ##

class predict_computational_cost:

    @staticmethod
    def ML_numevals(N_list):
        res_list = []
        for n in range(len(N_list) - 1):
            res_list.append(N_list[n] + N_list[n+1])
        res_list.append(N_list[-1])
        return res_list


    @staticmethod
    def NML_numevals(M_array):
        res_array = np.zeros(M_array.shape)
        for i in range(M_array.shape[0] - 1):
            res_array[i, 0] = np.sum(M_array[i:(i+2),:])
            res_array[i, 1] = np.sum(M_array[i:(i+2),1])
        res_array[-1, 0] = np.sum(M_array[-1,:])
        res_array[-1, 1] = M_array[-1, 1]
        return res_array


    @staticmethod
    def MLMF_numevals(N_list, r_list):
        res_list_hi = []
        res_list_lo = []

        for n in range(len(N_list) - 1):
            res_list_hi.append(N_list[n] + N_list[n+1])
            res_list_lo.append(N_list[n] + N_list[n+1] + r_list[n])

        res_list_hi.append(N_list[-1])
        res_list_lo.append(N_list[-1] + r_list[-1])

        return res_list_hi, res_list_lo


    @staticmethod
    def ML_time(numeval_list, runtime_list):
        total = 0
        for i in range(len(numeval_list)):
            total += numeval_list[i] * runtime_list[i]
        return total


    @staticmethod
    def NML_time(numeval_array, runtime_array):
        return np.sum(numeval_array * runtime_array)


    @staticmethod
    def MLMF_time(numeval_list_hi, numeval_list_lo, runtime_array):
        return ML_time(numeval_list_hi, runtime_array[:,0]) + ML_time(numeval_list_lo, runtime_array[:,1])
    
##################################################################################################

## INTERPOLATION ##

class interpolate:
    """Contains methods to interpolate SWAN output data on regular grids.
    
    - interpolate_time:     interpolate selected output variables to a finer time range;
    - get_nan_mask:         generates a mask of np.nan-values for two-dimensional value maps;
    - smoothly_filled_map:  fill the np.nan-valued region of a two-dimensional value map using biharmonic inpainting;
    - interpolate_space:    interpolate selected output variables to a finer spatial resolution, taking into account 
                            different np.nan-valued regions.
    """

    @staticmethod
    def interpolate_time(t, fp, output_vars=None):
        """Given an xArray dataset (e.g. SWAN output), interpolates all data variables temporally. Time ranges are assumed to contain np.datetime objects.
        
        - t (np.ndArray of np.datetime64 objects):      time range the data should be interpolated to.
        - fp (xArray dataset):                          data that should be interpolated temporally. Assumes dimensions to be 'time', 'longitude' and 'latitude'.
        - output_vars (list):                           list of output variables one wishes to interpolate. If output_vars is None, then all variables will be interpolated.
        """
        new_time = xr.DataArray(t, dims=('time'), name='time')
        new_lon = xr.DataArray(fp['longitude'].values, dims=('longitude'), name='longitude')
        new_lat = xr.DataArray(fp['latitude'].values, dims=('latitude'), name='latitude')

        new_time.attrs = fp['time'].attrs
        new_lon.attrs = fp['longitude'].attrs
        new_lat.attrs = fp['latitude'].attrs

        # Generate numerical timeranges

        num_new_time = np.array(t.values.tolist()) / 1000000000
        num_old_time = np.array(fp['time'].values.tolist()) / 1000000000

        new_time_grid, new_lat_grid, new_lon_grid = np.meshgrid(num_new_time, new_lat.values, new_lon.values, indexing='ij')

        new_datavars = {}

        if output_vars is None:
            output_vars = fp.keys()

        # Interpolate output variables
        
        for variable in output_vars:
            interp = RegularGridInterpolator((num_old_time, new_lat.values, new_lon.values), fp[variable].values[:,:,:])
            vals = interp((new_time_grid, new_lat_grid, new_lon_grid))      

            new_var = xr.DataArray(vals, dims=('time', 'latitude', 'longitude'), coords={'time': new_time, 'latitude': new_lat, 'longitude': new_lon}, name=variable)
            new_datavars[variable] = new_var

        interpolated_dataset = xr.Dataset(copy.deepcopy(new_datavars))
        interpolated_dataset.attrs = fp.attrs
        return interpolated_dataset


    @staticmethod
    def get_nan_mask(valmap):
        """Returns numpy array indicating with ones where np.nan-values are present.
        
        - valmap (np.ndarray):  array containing data with np.nan-values.
        """
        return np.where(~np.isnan(valmap), np.zeros(valmap.shape), 1)


    @staticmethod
    def smoothly_filled_map(valmap, nan_mask):
        """Fills np.nan-regions of a value map using biharmonic inpainting (scikit-image package).
        
        - valmap (np.ndarray):      array containing data with np.nan-values;
        - nan_mask (np.ndarray):    array indicating with ones, the location of np.nan-values (used as input of the function for computational speed).
        
        """
        filled_valmap = np.where(~np.isnan(valmap), valmap, 0)
        smoothed_valmap = inpaint.inpaint_biharmonic(filled_valmap, nan_mask)
        return smoothed_valmap
        

    @staticmethod
    def interpolate_space(new_fp, fp, output_vars=None):
        """Given an xArray dataset with np.nan-values (e.g. SWAN output), interpolates selected output variables spatially.
        
        - new_fp (xArray Dataset):  dataset containing higher resolution longitude and latitude data. This dataset is also used to obtain nan-masks
                                    for output-variables.
        - fp (xArray Dataset):      dataset containing the to be interpolated data.
        - output_vars (list):       list of output variables that will be interpolated. If output_vars is None, then all variables will be interpolated. 
        """
        new_time = xr.DataArray(fp['time'].values, dims=('time'), name='time')
        new_lon_data = xr.DataArray(new_fp['longitude'].values, dims=('longitude'), name='longitude')
        new_lat_data = xr.DataArray(new_fp['latitude'].values, dims=('latitude'), name='latitude')

        new_time.attrs = fp['time'].attrs
        new_lon_data.attrs = fp['longitude'].attrs
        new_lat_data.attrs = fp['latitude'].attrs

        new_lon_grid, new_lat_grid = np.meshgrid(new_fp['longitude'].values, new_fp['latitude'].values)

        new_datavars = {}

        if output_vars is None:
            output_vars = new_fp.keys()

        for variable in output_vars:
            new_data = np.zeros((new_time.values.shape[0], new_lat_data.values.shape[0], new_lon_data.values.shape[0]))
            nan_mask = interpolate.get_nan_mask(fp[variable][0,:,:].values)
            new_nan_mask = interpolate.get_nan_mask(new_fp[variable][0,:,:].values)
            for t in range(new_time.values.shape[0]):
                inpainted_map = interpolate.smoothly_filled_map(fp[variable][t,:,:].values, nan_mask) #inpaint the original map so edges along the land boundary blur
                
                interp = RegularGridInterpolator((fp['latitude'].values, fp['longitude'].values), inpainted_map) # interpolate the resulting array
                new_map = interp((new_lat_grid, new_lon_grid))

                new_data[t,:,:] = np.where(new_nan_mask == 0, new_map, np.nan) # set the correct coordinates to nan

            
            new_datavars[variable] = xr.DataArray(new_data, dims=('time', 'latitude', 'longitude'), coords={'time': new_time, 'longitude': new_lon_data, 'latitude': new_lat_data}, name=variable)
            new_datavars[variable].attrs = fp[variable].attrs
        
        interpolated_dataset = xr.Dataset(copy.deepcopy(new_datavars))
        interpolated_dataset.attrs = fp.attrs
        return interpolated_dataset

    
## MULTILEVEL/MULTIFIDELITY ESTIMATORS ##

class estimate:
    """Contains methods to estimate statistics of multilevel ensembles of SWAN output on regular grids.

    - ensemble_mean_dataset:            constructs an xArray dataset of the ensemble mean of an iid ensemble;   
    - ML_mean_estimator:                estimate the ensemble mean using a multilevel ensemble with the multilevel Monte Carlo method (MLMC);
    - variance_estimator:               estimate the ensemble variance using a multilevel ensemble with any mean approximation method;
    - ML_representative_ensemble:       using the method by Gregory and Cotter (2017), generates an approximate high-level iid ensemble at one
                                        location and one point in time; from this, quantile time series or maps may be computed;
    - quantile_from_ensemble:           estimate a quantile from an ensemble at a single point (full or approximated);
    - empirical_cdf:                    compute an empirical distribution function from an ensemble at a single point (full or approximated);
    - used_sample_indices_MLMF:         for a given approximation level index, determine which input sample indices have been used for simulations
                                        of this approximation level in the multilevel multifidelity Monte Carlo method (MLMF);
    - get_low_fidelity_expectations:    computes expected values of low-fidelity model approximations for the control variate formulation of MLMF;
    - MLMF_mean_estimator:              estimate the ensemble mean using a multilevel ensemble with the multilevel multifidelity Monte Carlo method (MLMF);
    - MLMF_representative_ensemble:     using the method by Clare et al. (2022), generates an approximate high-level high-fidelity iid ensemble
                                        at one point in space and time.   
    - NML_mean_estimator:               estimate the ensemble mean using a multilevel ensemble with the nested multilevel Monte Carlo method;
    - NML_representative_ensemble:      using a method analogous to Gregory & Cotter (2017), generates an approximate high-level iid ensemble at
                                        one point in space and time. 
    """

    @staticmethod
    def ensemble_mean_dataset(path, name, ensemble_size, output_vars = None, squared = False):
        """Returns xArray dataset of the ensemble mean of an iid ensemble. Can also be other data than SWAN data, such as TIGGE-wind.
        
        - path (str):           path to folder containing all ensemble members;
        - name (str):           name convention of ensemble members (without member index and extension);
        - output_vars (list):   list of output_variables that should be computed;
        - squared (bool):       if True, then the mean of the square will also be computed (for purposes of variance calculation).
        """
        # Open reference dataset
        first_member = xr.open_dataset(path + name + '0.nc', engine='netcdf4')

        if output_vars is None:
            output_vars = first_member.keys()

        # Initialize ensemble mean data

        data_coords = dict()
        for dim in first_member.dims:
            data_coords[dim] = first_member[dim]

        data = {}
        for variable in output_vars:
            data[variable] = xr.DataArray(np.zeros(first_member[variable].values.shape), dims = tuple(first_member.dims), coords = data_coords, name=variable)
            data[variable].attrs = first_member[variable].attrs

        ensemble_mean = xr.Dataset(copy.deepcopy(data))
        ensemble_mean.attrs = first_member.attrs

        if squared:
            ensemble_mean_of_square = xr.Dataset(copy.deepcopy(data))
            ensemble_mean_of_square.attrs = first_member.attrs

        # Add first member data to ensemble mean

        for variable in output_vars:
            ensemble_mean[variable].values += first_member[variable].values / ensemble_size
            if squared:
                ensemble_mean_of_square[variable].values += np.power(first_member[variable].values, 2) / ensemble_size

        # Loop over the rest of the ensemble

        for n in range(1, ensemble_size):
            member = xr.open_dataset(path + name + f'{n}.nc', engine='netcdf4')
            for variable in output_vars:
                ensemble_mean[variable].values += member[variable].values / ensemble_size
                if squared:
                    ensemble_mean_of_square[variable].values += np.power(member[variable].values, 2) / ensemble_size

        # Return
        
        if squared:
            return ensemble_mean, ensemble_mean_of_square
        else:
            return ensemble_mean



    @staticmethod
    def ML_mean_estimator(path, run_name, N_list, reference_data, time_interpolation, space_interpolation, output_vars = None, squared = False):
        """Returns xArray dataset of the multilevel-estimated ensemble mean, where the data variables are given in output_vars.
        
        - path (str):                       path to folder containing all the EPS runs,
        - run_name (str):                   name of this EPS run;
        - N_list (list):                    values of N used for the MLMC method.
        - time_interpolation (list):        list of booleans indicating for each level whether they have to be interpolated in time;
        - space_interpolation (list):       list of booleans indicating for each level whether they have to be interpolated in time;
        - reference_data (xArray dataset):  dataset to reference for interpolation and metadata purposes. Should be a highest-level ensemble member;
        - output_vars (list):               list of variables to estimate the ensemble mean for;
        - squared (bool):                   if True, then the mean of the square will also be computed (for purposes of variance calculation).
        """

        # initialisation

        time = xr.DataArray(reference_data['time'].values, dims=('time'), name='time')
        latitude = xr.DataArray(reference_data['latitude'].values, dims=('latitude'), name='latitude')
        longitude = xr.DataArray(reference_data['longitude'].values, dims=('longitude'), name='longitude')

        time.attrs = reference_data['time'].attrs
        latitude.attrs = reference_data['latitude'].attrs
        longitude.attrs = reference_data['longitude'].attrs

        if output_vars is None:
            output_vars = reference_data.keys()

        data = {} 
        for variable in output_vars:
            data[variable] = xr.DataArray(np.zeros(reference_data[variable].values.shape), dims=('time', 'latitude', 'longitude'), coords={'time':time, 'latitude': latitude, 'longitude': longitude}, name=variable)
            data[variable].attrs = reference_data[variable].attrs

        ensemble_mean = xr.Dataset(copy.deepcopy(data))
        ensemble_mean.attrs = reference_data.attrs

        if squared:
            ensemble_mean_of_square = xr.Dataset(copy.deepcopy(data))
            ensemble_mean_of_square.attrs = reference_data.attrs


        # LEVEL 0 #

        for i in range(N_list[0]):
            sample_number = i
            # load ensemble member
            member = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level0_sample{sample_number}.nc', engine='netcdf4')

            # interpolate if necessary
            if time_interpolation[0]:
                t_interpolated_member = interpolate.interpolate_time(reference_data['time'].values, member, output_vars=output_vars)
            else:
                t_interpolated_member = member

            if space_interpolation[0]:
                interpolated_member = interpolate.interpolate_space(reference_data, t_interpolated_member, output_vars=output_vars)
            else:
                interpolated_member = t_interpolated_member

            # add to dataset
            
            for variable in output_vars:
                ensemble_mean[variable].values += interpolated_member[variable].values / N_list[0]
                if squared:
                    ensemble_mean_of_square[variable].values += np.power(interpolated_member[variable].values, 2) / N_list[0]


        # HIGHER LEVEL CORRECTIONS #

        for l in range(1, len(N_list)):
            for i in range(N_list[l]):
                sample_number = sum(N_list[:l]) + i
                member_hi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_sample{sample_number}.nc', engine='netcdf4')
                member_lo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l-1}_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[l]:
                    t_interpolated_member_hi = interpolate.interpolate_time(reference_data['time'].values, member_hi, output_vars=output_vars)
                else:
                    t_interpolated_member_hi = member_hi

                if time_interpolation[l-1]:
                    t_interpolated_member_lo = interpolate.interpolate_time(reference_data['time'].values, member_lo, output_vars=output_vars)
                else:
                    t_interpolated_member_lo = member_lo

                if space_interpolation[l]:
                    interpolated_member_hi = interpolate.interpolate_space(reference_data, t_interpolated_member_hi, output_vars=output_vars)
                else:
                    interpolated_member_hi = t_interpolated_member_hi

                if space_interpolation[l-1]:
                    interpolated_member_lo = interpolate.interpolate_space(reference_data, t_interpolated_member_lo, output_vars=output_vars)
                else:
                    interpolated_member_lo = t_interpolated_member_lo         

                # add to ensemble mean dataset
                for variable in output_vars:
                    ensemble_mean[variable].values += (interpolated_member_hi[variable].values - interpolated_member_lo[variable].values) / N_list[l]
                    if squared:
                        ensemble_mean_of_square[variable].values += (np.power(interpolated_member_hi[variable].values, 2) - \
                                                                np.power(interpolated_member_lo[variable].values, 2)) / N_list[l]
        if squared:
            return ensemble_mean, ensemble_mean_of_square
        else:
            return ensemble_mean


    @staticmethod
    def variance_estimator(ensemble_mean, ensemble_mean_of_square, output_vars=None):
        """Returns xArray dataset of the multilevel-estimated ensemble variance (any method), where the data variables are given in output_vars.

        - ensemble_mean (xArray dataset):           precomputed ensemble mean dataset;
        - ensemble_mean_of_square (xArray dataset): precomputed ensemble mean dataset of the square ensemble (output of *_mean_estimator with
                                                    squared==True);
        - output_vars (list):                       list of variables to estimate the ensemble mean for.
        """
        
        data = {}
        
        if output_vars is None:
            output_vars = ensemble_mean.keys()

        # Initialise dataset

        data_coords = dict()
        for dim in ensemble_mean.dims:
            data_coords[dim] = ensemble_mean[dim].values

        for variable in output_vars:
            data[variable] = xr.DataArray(np.zeros(ensemble_mean[variable].values.shape), dims=ensemble_mean.dims,
                                     coords = data_coords)
        
        ensemble_variance = xr.Dataset(copy.deepcopy(data))
        ensemble_variance.attrs = ensemble_mean.attrs

        # Compute the variance
        for variable in output_vars:
            ensemble_variance[variable].values = ensemble_mean_of_square[variable].values - np.power(ensemble_mean[variable].values, 2)

        return ensemble_variance
    

    @staticmethod
    def ML_representative_ensemble(path, run_name, N_list, reference_data, time_interpolation, space_interpolation, longitude_point, latitude_point, time_index, output_var, size):
        """Generates approximate high-level ensemble at one point in space and time with the inverse transform sampling method described in
        Gregory and Cotter (2017). From this representative ensemble, empirical distribution functions and quantiles may be computed.
        
        - path (str):                       path to folder containing all the EPS runs;
        - run_name (str):                   name of the EPS run of interest;
        - N_list (list):                    list of numbers N used in this MLMC EPS run;
        - reference_data (xArray dataset):  dataset used to know how interpolation should be done;
        - time_interpolation (list):        list of booleans indicating whether each level should be interpolated in time;
        - space_interpolation (list):       list of booleans indicating whether each level should be interpolated in space;
        - longitude_point (float):          longitude of the point of interest;
        - latitude_point (float):           latitude of the point of interest;
        - time_index (int):                 time index of interest;
        - output_var (str):                 output variable of interest;
        - size (int):                       size of the to be generated ensemble.
        """
        # Compute latitude and longitude indices
        lat_index = np.argmin(np.absolute(np.full(reference_data['latitude'].values.shape, latitude_point) - reference_data['latitude'].values))
        lon_index = np.argmin(np.absolute(np.full(reference_data['longitude'].values.shape, longitude_point) - reference_data['longitude'].values))
        point_of_interest = (time_index, lat_index, lon_index)

        # Gather data at point of interest into a numpy array
        point_ensemble = np.zeros(N_list[0] + 2 * sum(N_list[1:]))

        for i in range(N_list[0]):
            sample_number = i
            member = xr.open_dataset(path+f'{run_name}\\output\\{run_name}_level0_sample{sample_number}.nc', engine='netcdf4')
            
            # interpolate if necessary
            if time_interpolation[0]:
                t_interpolated_member = interpolate.interpolate_time(reference_data['time'].values, member, output_vars=[output_var])
            else:
                t_interpolated_member = member

            if space_interpolation[0]:
                interpolated_member = interpolate.interpolate_space(reference_data, t_interpolated_member, output_vars=[output_var])
            else:
                interpolated_member = t_interpolated_member

            point_ensemble[i] = interpolated_member[output_var].values[point_of_interest]

        for l in range(1, len(N_list)):
            for i in range(N_list[l]):
                sample_number = sum(N_list[:l]) + i
                member_hi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_sample{sample_number}.nc', engine='netcdf4')
                member_lo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l-1}_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[l]:
                    t_interpolated_member_hi = interpolate.interpolate_time(reference_data['time'].values, member_hi, output_vars=[output_var])
                else:
                    t_interpolated_member_hi = member_hi

                if time_interpolation[l-1]:
                    t_interpolated_member_lo = interpolate.interpolate_time(reference_data['time'].values, member_lo, output_vars=[output_var])
                else:
                    t_interpolated_member_lo = member_lo

                if space_interpolation[l]:
                    interpolated_member_hi = interpolate.interpolate_space(reference_data, t_interpolated_member_hi, output_vars=[output_var])
                else:
                    interpolated_member_hi = t_interpolated_member_hi

                if space_interpolation[l-1]:
                    interpolated_member_lo = interpolate.interpolate_space(reference_data, t_interpolated_member_lo, output_vars=[output_var])
                else:
                    interpolated_member_lo = t_interpolated_member_lo   

                point_ensemble[N_list[0] + 2*sum(N_list[1:l]) + i] = interpolated_member_lo[output_var].values[point_of_interest]
                point_ensemble[N_list[0] + 2*sum(N_list[1:l]) + N_list[l] + i] = interpolated_member_hi[output_var].values[point_of_interest]

        # Sort the ensemble in the correct way (described in report)

        sorted_ensemble = np.zeros(N_list[0] + 2*sum(N_list[1:]))

        sorted_ensemble[:N_list[0]] = np.sort(point_ensemble[:N_list[0]])

        for l in range(1, len(N_list)):
            start_index = N_list[0] + 2*sum(N_list[1:l])
            mid_index = start_index + N_list[l]
            end_index = mid_index + N_list[l]

            sorted_ensemble[start_index:mid_index] = np.sort(point_ensemble[start_index:mid_index])
            sorted_ensemble[mid_index:end_index] = np.sort(point_ensemble[mid_index:end_index])

        # Generate representative ensemble using inverse transform sampling method by Gregory & Cotter (2017)

        repr_ensemble = np.zeros(size)

        for n in range(size):
            u = np.random.uniform()
            repr_ensemble[n] += sorted_ensemble[int(np.floor(N_list[0]*u))] # use floor instead of ceil (in the paper), since we count from zero

            for l in range(1, len(N_list)):
                start_index = N_list[0] + 2*sum(N_list[1:l])
                mid_index = start_index + N_list[l]

                repr_ensemble[n] += sorted_ensemble[mid_index + int(np.floor(N_list[l]*u))] - \
                                    sorted_ensemble[start_index + int(np.floor(N_list[l]*u))]
    
        return repr_ensemble


    @staticmethod
    def quantile_from_ensemble(u, ensemble):
        """Estimates the u-quantile of a scalar probability distribution using a scalar ensemble.
        
        - u (float):                quantile (e.g. u=0.5 gives the median);
        - ensemble (np.ndarray):    scalar ensemble assumed to be independent and identically distributed.
        """

        sorted_ensemble = np.sort(ensemble)
        ensemble_size = ensemble.shape[0]
        return sorted_ensemble[int(np.floor(ensemble_size * u))]
    

    @staticmethod
    def empirical_cdf(num_steps, ensemble):
        """Compute the empirical distribution function based on a scalar ensemble. Returns coordinate array x (x-axis of the cdf) and corresponding
        distribution function values.
        
        - num_steps (int):          step size of the x axis;
        - ensemble (np.ndarray):    scalar ensemble assumed to be independent and identically distributed.
        """
        size = ensemble.shape[0]
        xaxis = np.linspace(np.amin(ensemble), np.amax(ensemble), num_steps)
        distribution_function = np.zeros(num_steps)

        for i in range(num_steps):
            # if X^(i) < x, then sign(x-X^(i))=1
            indicator = np.maximum(np.sign(np.full(size, xaxis[i]) - ensemble), np.zeros(size))
            distribution_function[i] = np.sum(indicator) / size

        return xaxis, distribution_function
    

    @staticmethod
    def used_sample_indices_MLMF(N_list, r_list, level_index):
        """Returns a list of indices specifying which input samples were used to simulate the collection of low-fidelity runs.

        - N_list (list):        values of N^(HF) used in the MLMF-method;
        - r_list (list):        amounts of extra simulations performed to estimate the mean of the low-fidelity models;
        - level_index (int):    index specifying which approximation level should be considered;
        """

        index_list = []

        # Samples also used for high-fidelity model

        if level_index < len(N_list) - 1:
            for i in range(N_list[level_index] +  N_list[level_index+1]):
                index_list.append(sum(N_list[:level_index]) + i)
        else:
            for i in range(N_list[level_index]):
                index_list.append(sum(N_list[:level_index]) + i)

        # Extra samples

        input_start_index = sum(N_list[:level_index])

        if level_index == 0:
            for i in range(r_list[level_index]):
                index_list.append(N_list[0] + N_list[1] + i)
        else:
            # Indices before already used ones
            for i in range(min(r_list[level_index], input_start_index)):
                index_list.append(i)

            # Indices after already used ones
            for i in range(max(0, r_list[level_index] - input_start_index)):
                if level_index < len(N_list) - 1:
                    index_list.append(input_start_index + N_list[level_index] + N_list[level_index + 1] + i)
                else:
                    index_list.append(input_start_index + N_list[level_index] + i)

        return index_list


    @staticmethod
    def get_low_fidelity_expectations(path, run_name, N_list, r_list, reference_data, time_interpolation, space_interpolation, output_vars = None, squared = False):
        """Returns a list of xArrays containing the means of the low-fidelity model at each approximation level.
        
        - path (str):                       path to folder containing all the EPS runs;
        - run_name (str):                   name of the EPS run of interest;
        - N_list (list):                    values of N^(HF) used in the MLMF-method;
        - r_list (list):                    amounts of extra simulations performed to estimate the mean of the low-fidelity models;
        - time_interpolation (np.ndarray):  array of booleans indicating whether temporal interpolation should take place. Column index 
                                            indicates approximation level and row index specifies high (0) or low (1) fidelity;
        - space_interpolation (np.ndarray): array of booleans indicating whether spatial interpolation should take place. Column index
                                            indicates approximation level and row index specifies high (0) or low (1) fidelity;
        - reference_data (xArray dataset):  dataset used to gather metadata and fine grids for interpolation purposes;
        - output_vars (list):               selection of output variables that will be computed;
        - squared (bool):                   if True, then the mean of the square ensemble will be estimated too. This is useful for 
                                            variance estimates.
        """

        # Initialise data

        time = xr.DataArray(reference_data['time'].values, dims=('time'), name='time')
        latitude = xr.DataArray(reference_data['latitude'].values, dims=('latitude'), name='latitude')
        longitude = xr.DataArray(reference_data['longitude'].values, dims=('longitude'), name='longitude')

        time.attrs = reference_data['time'].attrs
        latitude.attrs = reference_data['latitude'].attrs
        longitude.attrs = reference_data['longitude'].attrs

        if output_vars is None:
            output_vars = reference_data.keys()

        data = {} 
        for variable in output_vars:
            data[variable] = xr.DataArray(np.zeros(reference_data[variable].values.shape), dims=('time', 'latitude', 'longitude'), 
                                     coords={'time':time, 'latitude': latitude, 'longitude': longitude}, name=variable)
            data[variable].attrs = reference_data[variable].attrs

        expectations = []
        for l in range(len(N_list)):
            expectations.append(xr.Dataset(copy.deepcopy(data)))
            expectations[l].attrs = reference_data.attrs

        if squared:
            expectations_of_square = []
            for l in range(len(N_list)):
                expectations_of_square.append(xr.Dataset(copy.deepcopy(data)))
                expectations_of_square[l].attrs = reference_data.attrs
        
        # Compute expectations

        for l in range(len(N_list)):
            index_list = estimate.used_sample_indices_MLMF(N_list, r_list, l)

            for i in index_list:
                member = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_lo_sample{i}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[l, 1]:
                    t_interpolated_member = interpolate.interpolate_time(reference_data['time'].values, member, output_vars=output_vars)
                else:
                    t_interpolated_member = member

                if space_interpolation[l, 1]:
                    interpolated_member = interpolate.interpolate_space(reference_data, t_interpolated_member, output_vars=output_vars)
                else:
                    interpolated_member = t_interpolated_member

                # add data to mean computation

                for variable in output_vars:
                    expectations[l][variable].values += interpolated_member[variable].values / len(index_list)

                if squared:
                    for variable in output_vars:
                        expectations_of_square[l][variable].values += np.power(interpolated_member[variable].values, 2) / len(index_list)
                

        if squared:
            return expectations, expectations_of_square
        else:
            return expectations


    @staticmethod
    def MLMF_mean_estimator(path, run_name, N_list, r_list, alpha_list, reference_data, time_interpolation, space_interpolation, output_vars = None, squared = False):
        """Returns xArray dataset of the estimated ensemble mean using the multilevel multifidelity Monte Carlo method (MLMF).
        
        - path (str):                       path to folder containing all the EPS runs;
        - run_name (str):                   name of the EPS run of interest;
        - N_list (list):                    values of N^(HF) used in the MLMF-method;
        - r_list (list):                    amounts of extra simulations performed to estimate the mean of the low-fidelity models;
        - alpha_list (list):                list of parameters alpha used in the MLMF-method;
        - time_interpolation (np.ndarray):  array of booleans indicating whether temporal interpolation should take place. Row index 
                                            indicates approximation level and column index specifies high (0) or low (1) fidelity;
        - space_interpolation (np.ndarray): array of booleans indicating whether spatial interpolation should take place. Row index
                                            indicates approximation level and column index specifies high (0) or low (1) fidelity;
        - reference_data (xArray dataset):  dataset used to gather metadata and fine grids for interpolation purposes;
        - output_vars (list):               selection of output variables that will be computed;
        - squared (bool):                   if True, then the mean of the square ensemble will be estimated too. This is useful for 
                                            variance estimates.
        """

        # Initialise data

        time = xr.DataArray(reference_data['time'].values, dims=('time'), name='time')
        latitude = xr.DataArray(reference_data['latitude'].values, dims=('latitude'), name='latitude')
        longitude = xr.DataArray(reference_data['longitude'].values, dims=('longitude'), name='longitude')

        time.attrs = reference_data['time'].attrs
        latitude.attrs = reference_data['latitude'].attrs
        longitude.attrs = reference_data['longitude'].attrs

        if output_vars is None:
            output_vars = reference_data.keys()

        data = {} 
        for variable in output_vars:
            data[variable] = xr.DataArray(np.zeros(reference_data[variable].values.shape), dims=('time', 'latitude', 'longitude'), coords={'time':time, 'latitude': latitude, 'longitude': longitude}, name=variable)
            data[variable].attrs = reference_data[variable].attrs

        ensemble_mean = xr.Dataset(copy.deepcopy(data))
        ensemble_mean.attrs = reference_data.attrs

        if squared:
            ensemble_mean_of_square = xr.Dataset(copy.deepcopy(data))
            ensemble_mean_of_square.attrs = reference_data.attrs

        # Get expectations

        if squared:
            expectations, expectations_of_square = estimate.get_low_fidelity_expectations(path, run_name, N_list, r_list, time_interpolation,
                                                                                          space_interpolation, reference_data, output_vars=output_vars,
                                                                                          squared=True)
        else:
            expectations = estimate.get_low_fidelity_expectations(path, run_name, N_list, r_list, time_interpolation, space_interpolation, reference_data,
                                                                  output_vars=output_vars, squared=False)
            
        
        ## LEVEL 0 ##

        for i in range(N_list[0]):
            sample_number = i
            member_hf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level0_hi_sample{sample_number}.nc', engine='netcdf4')
            member_lf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level0_lo_sample{sample_number}.nc', engine='netcdf4')

            # interpolate if necessary

            if time_interpolation[0, 0]:
                t_interpolated_member_hf = interpolate.interpolate_time(reference_data['time'].values, member_hf, output_vars=output_vars)
            else:
                t_interpolated_member_hf = member_hf

            if time_interpolation[0, 1]:
                t_interpolated_member_lf = interpolate.interpolate_time(reference_data['time'].values, member_lf, output_vars=output_vars)
            else:
                t_interpolated_member_lf = member_lf

            if space_interpolation[0, 0]:
                interpolated_member_hf = interpolate.interpolate_space(reference_data, t_interpolated_member_hf, output_vars=output_vars)
            else:
                interpolated_member_hf = t_interpolated_member_hf
            
            if space_interpolation[0, 1]:
                interpolated_member_lf = interpolate.interpolate_space(reference_data, t_interpolated_member_lf, output_vars=output_vars)
            else:
                interpolated_member_lf = t_interpolated_member_lf

            # add data to ensemble mean computation

            for variable in output_vars:
                # it is assumed in this method that the parameter alpha is the same for all output variables
                ensemble_mean[variable].values += (interpolated_member_hf[variable].values + alpha_list[0] * interpolated_member_lf[variable].values) / N_list[0]

                if squared:
                    ensemble_mean_of_square[variable].values += (np.power(interpolated_member_hf[variable].values, 2) + alpha_list[0] * \
                                                            np.power(interpolated_member_lf[variable].values, 2)) / N_list[0]
            
        
        ## HIGHER LEVEL CORRECTIONS ##

        for l in range(1, len(N_list)):
            for i in range(N_list[l]):
                sample_number = sum(N_list[:l]) + i

                # load data
                member_hlhf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_hi_sample{sample_number}.nc', engine='netcdf4')
                member_hllf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_lo_sample{sample_number}.nc', engine='netcdf4')
                member_llhf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l-1}_hi_sample{sample_number}.nc', engine='netcdf4')
                member_lllf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l-1}_lo_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[l, 0]:
                    t_interpolated_member_hlhf = interpolate.interpolate_time(reference_data['time'].values, member_hlhf, output_vars=output_vars)
                else:
                    t_interpolated_member_hlhf = member_hlhf

                if time_interpolation[l, 1]:
                    t_interpolated_member_hllf = interpolate.interpolate_time(reference_data['time'].values, member_hllf, output_vars=output_vars)
                else:
                    t_interpolated_member_hllf = member_hllf

                if space_interpolation[l, 0]:
                    interpolated_member_hlhf = interpolate.interpolate_space(reference_data, t_interpolated_member_hlhf, output_vars=output_vars)
                else:
                    interpolated_member_hlhf = t_interpolated_member_hlhf
                
                if space_interpolation[l, 1]:
                    interpolated_member_hllf = interpolate.interpolate_space(reference_data, t_interpolated_member_hllf, output_vars=output_vars)
                else:
                    interpolated_member_hllf = t_interpolated_member_hllf

                if time_interpolation[l-1, 0]:
                    t_interpolated_member_llhf = interpolate.interpolate_time(reference_data['time'].values, member_llhf, output_vars=output_vars)
                else:
                    t_interpolated_member_llhf = member_llhf

                if time_interpolation[l-1, 1]:
                    t_interpolated_member_lllf = interpolate.interpolate_time(reference_data['time'].values, member_lllf, output_vars=output_vars)
                else:
                    t_interpolated_member_lllf = member_lllf

                if space_interpolation[l-1, 0]:
                    interpolated_member_llhf = interpolate.interpolate_space(reference_data, t_interpolated_member_llhf, output_vars=output_vars)
                else:
                    interpolated_member_llhf = t_interpolated_member_llhf
                
                if space_interpolation[l-1, 1]:
                    interpolated_member_lllf = interpolate.interpolate_space(reference_data, t_interpolated_member_lllf, output_vars=output_vars)
                else:
                    interpolated_member_lllf = t_interpolated_member_lllf

                # add data to ensemble mean computation

                for variable in output_vars:
                    ensemble_mean[variable].values += (1 / N_list[l]) * (
                        (interpolated_member_hlhf[variable].values - interpolated_member_llhf[variable].values) + \
                        alpha_list[l] * (interpolated_member_hllf[variable].values - interpolated_member_lllf[variable].values - \
                                         expectations[l][variable].values + expectations[l-1][variable].values)
                    )

                    if squared:
                        ensemble_mean_of_square[variable].values += (1 / N_list[l]) * (
                            np.power(interpolated_member_hlhf[variable].values, 2) - np.power(interpolated_member_llhf[variable].values, 2) + \
                            alpha_list[l] * (np.power(interpolated_member_hllf[variable].values, 2) - np.power(interpolated_member_lllf[variable].values, 2) - \
                                            expectations_of_square[l][variable].values + expectations_of_square[l-1][variable].values)
                        )
            
            if squared:
                return ensemble_mean, ensemble_mean_of_square
            else:
                return ensemble_mean
        
        
    @staticmethod
    def MLMF_representative_ensemble(path, run_name, N_list, expectations, alpha_list, reference_data, time_interpolation, space_interpolation, longitude_point, latitude_point, time_index, output_var, size):
        """Generates approximate high-level ensemble at one point in space and time with the inverse transform sampling method described in
        Clare et al. (2022). From this representative ensemble, empirical distribution functions and quantiles may be computed.
        
        - path (str):                       path to folder containing all the EPS runs;
        - run_name (str):                   name of the EPS run of interest;
        - N_list (list):                    list of numbers N used in this MLMC EPS run;
        - expectations (xArray dataset):    list of expected values of the low-fidelity model at every approximation level;
        - alpha_list (list):                list of parameters alpha used in MLMF-method;
        - reference_data (xArray dataset):  dataset used to know how interpolation should be done;
        - time_interpolation (list):        list of booleans indicating whether each level should be interpolated in time;
        - space_interpolation (list):       list of booleans indicating whether each level should be interpolated in space;
        - longitude_point (float):          longitude of the point of interest;
        - latitude_point (float):           latitude of the point of interest;
        - time_index (int):                 time index of interest;
        - output_var (str):                 output variable of interest;
        - size (int):                       size of the to be generated ensemble.
        """

        # Compute latitude and longitude indices
        lat_index = np.argmin(np.absolute(np.full(reference_data['latitude'].values.shape, latitude_point) - reference_data['latitude'].values))
        lon_index = np.argmin(np.absolute(np.full(reference_data['longitude'].values.shape, longitude_point) - reference_data['longitude'].values))
        point_of_interest = (time_index, lat_index, lon_index)

        # Gather data at point of interest
        point_ensemble_hf = np.zeros(N_list[0] + 2 * sum(N_list[1:]))
        point_ensemble_lf = np.zeros(N_list[0] + 2 * sum(N_list[1:]))

        point_expectations = [expectations[l][output_var].values[point_of_interest] for l in range(len(N_list))]

        for i in range(N_list[0]):
            sample_number = i
            member_hf = xr.open_dataset(path+f'{run_name}\\output\\{run_name}_level0_hi_sample{sample_number}.nc', engine='netcdf4')
            member_lf = xr.open_dataset(path+f'{run_name}\\output\\{run_name}_level0_lo_sample{sample_number}.nc', engine='netcdf4')
            
            # interpolate if necessary
            if time_interpolation[0, 0]:
                t_interpolated_member_hf = interpolate.interpolate_time(reference_data['time'].values, member_hf, output_vars=[output_var])
            else:
                t_interpolated_member_hf = member_hf

            if space_interpolation[0, 0]:
                interpolated_member_hf = interpolate.interpolate_space(reference_data, t_interpolated_member_hf, output_vars=[output_var])
            else:
                interpolated_member_hf = t_interpolated_member_hf  

            if time_interpolation[0, 1]:
                t_interpolated_member_lf = interpolate.interpolate_time(reference_data['time'].values, member_lf, output_vars=[output_var])
            else:
                t_interpolated_member_lf = member_hf

            if space_interpolation[0, 1]:
                interpolated_member_lf = interpolate.interpolate_space(reference_data, t_interpolated_member_lf, output_vars=[output_var])
            else:
                interpolated_member_lf = t_interpolated_member_lf

            point_ensemble_hf[i] = interpolated_member_hf[output_var].values[point_of_interest]
            point_ensemble_lf[i] = interpolated_member_lf[output_var].values[point_of_interest]

        for l in range(1, len(N_list)):
            for i in range(N_list[l]):
                sample_number = sum(N_list[:l]) + i
                member_hlhf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_hi_sample{sample_number}.nc', engine='netcdf4')
                member_llhf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l-1}_hi_sample{sample_number}.nc', engine='netcdf4')
                member_hllf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l}_lo_sample{sample_number}.nc', engine='netcdf4')
                member_lllf = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_level{l-1}_lo_sample{sample_number}.nc', engine='netcdf4')


                # interpolate if necessary
                if time_interpolation[l, 0]:
                    t_interpolated_member_hlhf = interpolate.interpolate_time(reference_data['time'].values, member_hlhf, output_vars=[output_var])
                else:
                    t_interpolated_member_hlhf = member_hlhf

                if time_interpolation[l-1, 0]:
                    t_interpolated_member_llhf = interpolate.interpolate_time(reference_data['time'].values, member_llhf, output_vars=[output_var])
                else:
                    t_interpolated_member_llhf = member_llhf

                if space_interpolation[l, 0]:
                    interpolated_member_hlhf = interpolate.interpolate_space(reference_data, t_interpolated_member_hlhf, output_vars=[output_var])
                else:
                    interpolated_member_hlhf = t_interpolated_member_hlhf

                if space_interpolation[l-1, 0]:
                    interpolated_member_llhf = interpolate.interpolate_space(reference_data, t_interpolated_member_llhf, output_vars=[output_var])
                else:
                    interpolated_member_llhf = t_interpolated_member_llhf   

                if time_interpolation[l, 1]:
                    t_interpolated_member_hllf = interpolate.interpolate_time(reference_data['time'].values, member_hllf, output_vars=[output_var])
                else:
                    t_interpolated_member_hllf = member_hllf

                if time_interpolation[l-1, 1]:
                    t_interpolated_member_lllf = interpolate.interpolate_time(reference_data['time'].values, member_lllf, output_vars=[output_var])
                else:
                    t_interpolated_member_lllf = member_lllf

                if space_interpolation[l, 1]:
                    interpolated_member_hllf = interpolate.interpolate_space(reference_data, t_interpolated_member_hllf, output_vars=[output_var])
                else:
                    interpolated_member_hllf = t_interpolated_member_hllf

                if space_interpolation[l-1, 1]:
                    interpolated_member_lllf = interpolate.interpolate_space(reference_data, t_interpolated_member_lllf, output_vars=[output_var])
                else:
                    interpolated_member_lllf = t_interpolated_member_lllf  

                point_ensemble_hf[N_list[0] + 2*sum(N_list[1:l]) + i] = interpolated_member_hlhf[output_var].values[point_of_interest]
                point_ensemble_hf[N_list[0] + 2*sum(N_list[1:l]) + N_list[l] + i] = interpolated_member_llhf[output_var].values[point_of_interest]
                point_ensemble_lf[N_list[0] + 2*sum(N_list[1:l]) + i] = interpolated_member_hllf[output_var].values[point_of_interest]
                point_ensemble_lf[N_list[0] + 2*sum(N_list[1:l]) + N_list[l] + i] = interpolated_member_lllf[output_var].values[point_of_interest]

        # Sort the point ensembles in the correct way

        sorted_ensemble_hf = np.zeros(point_ensemble_hf.shape[0])
        sorted_ensemble_lf = np.zeros(point_ensemble_lf.shape[0])

        # Level 0
        sorted_ensemble_hf[:N_list[0]] = np.sort(point_ensemble_hf[:N_list[0]])
        sorted_ensemble_lf[:N_list[0]] = np.sort(point_ensemble_lf[:N_list[0]])

        # Higher levels
        for l in range(1, len(N_list)):
            start_index = N_list[0] + 2*sum(N_list[1:l])
            mid_index = start_index + N_list[l]
            end_index = mid_index + N_list[l]

            sorted_ensemble_hf[start_index:mid_index] = np.sort(point_ensemble_hf[start_index:mid_index])
            sorted_ensemble_hf[mid_index:end_index] = np.sort(point_ensemble_hf[mid_index:end_index])

            sorted_ensemble_lf[start_index:mid_index] = np.sort(point_ensemble_lf[start_index:mid_index])
            sorted_ensemble_lf[mid_index:end_index] = np.sort(point_ensemble_lf[mid_index:end_index])

        # Sample using the inverse transform sampling technique described in Clare et al. (2022)

        repr_ensemble = np.zeros(size)

        for n in range(size):
            u = np.random.uniform()

            repr_ensemble[n] = sorted_ensemble_hf[int(np.floor(N_list[0] * u))] # low-level high-fidelity

            repr_ensemble[n] += alpha_list[0] * (sorted_ensemble_lf[int(np.floor(N_list[0] * u))] - point_expectations[0]) # low-level low-fidelity correction
        
            # Higher level corrections
            
            for l in range(1, len(N_list)):
                start_index = N_list[0] + 2*sum(N_list[1:l])
                mid_index = start_index + N_list[l]

                repr_ensemble[n] += sorted_ensemble_hf[mid_index + int(np.floor(N_list[l]*u))] - \
                                        sorted_ensemble_hf[start_index + int(np.floor(N_list[l]*u))]
                
                repr_ensemble[n] += alpha_list[l] * (sorted_ensemble_lf[mid_index + int(np.floor(N_list[l]*u))] - sorted_ensemble_lf[start_index + int(np.floor(N_list[l]*u))]
                                                    - (point_expectations[l] - point_expectations[l-1]))                
            
        return repr_ensemble  
            

    @staticmethod
    def NML_mean_estimator(path, run_name, M_array, reference_data, time_interpolation, space_interpolation, output_vars = None, squared = False):
        """Returns xArray dataset of estimated ensemble mean using a multilevel ensemble with the nested multilevel Monte Carlo method (NML).
        
        - path (str):                           path to folder containing all the EPS runs;
        - run_name (str):                       name of the EPS run of interest;
        - M_array (str):                        values of M used in the NML method. Row index denotes outer approximation level and column index
                                                denotes inner approximation level;
        - reference_data (xArray dataset):      data used to obtain metadata and fine grids for interpolation purposes;
        - time_interpolation (np.ndarray):      array of booleans indicating whether model approximations should be interpolated in time. Row
                                                index denotes outer approximation level and column index denotes inner approximation level;
        - space_interpolation (np.ndarray):     array of booleans indicating whether model approximations should be interpolated in space. Row
                                                index denotes outer approximation level and column index denotes inner approximation level;
        - output_vars (list):                   list of output variables to estimate the mean of. If it is None, then all output variables will be used;
        - squared (bool):                       flag indicating whether the mean of the squared ensemble should be computed as well.
        """

        # Initialise data

        time = xr.DataArray(reference_data['time'].values, dims=('time'), name='time')
        latitude = xr.DataArray(reference_data['latitude'].values, dims=('latitude'), name='latitude')
        longitude = xr.DataArray(reference_data['longitude'].values, dims=('longitude'), name='longitude')

        time.attrs = reference_data['time'].attrs
        latitude.attrs = reference_data['latitude'].attrs
        longitude.attrs = reference_data['longitude'].attrs

        if output_vars is None:
            output_vars = reference_data.keys()

        data = {} 
        for variable in output_vars:
            data[variable] = xr.DataArray(np.zeros(reference_data[variable].values.shape), dims=('time', 'latitude', 'longitude'), coords={'time':time, 'latitude': latitude, 'longitude': longitude}, name=variable)
            data[variable].attrs = reference_data[variable].attrs

        ensemble_mean = xr.Dataset(copy.deepcopy(data))
        ensemble_mean.attrs = reference_data.attrs

        if squared:
            ensemble_mean_of_square = xr.Dataset(copy.deepcopy(data))
            ensemble_mean_of_square.attrs = reference_data.attrs

        num_outer_levels = M_array.shape[0]
        num_inner_levels = M_array.shape[1]

        ## OUTER LEVEL 0 ##

            # Inner level 0
        for i in range(M_array[0, 0]):
            sample_number = i
            member = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel0_ilevel0_sample{sample_number}.nc', engine='netcdf4')

            # interpolate if necessary
            if time_interpolation[0, 0]:
                t_interpolated_member = interpolate.interpolate_time(reference_data['time'].values, member, output_vars=output_vars)
            else:
                t_interpolated_member = member

            if space_interpolation[0, 0]:
                interpolated_member = interpolate.interpolate_space(reference_data, t_interpolated_member, output_vars=output_vars)
            else:
                interpolated_member = t_interpolated_member

            # add to ensemble mean

            for variable in output_vars:
                ensemble_mean[variable].values += interpolated_member[variable].values / M_array[0, 0]
                if squared:
                    ensemble_mean_of_square[variable].values += np.power(interpolated_member[variable].values, 2) / M_array[0,0]

            
            # Higher inner level corrections
        for k in range(1, num_inner_levels):
            for i in range(M_array[0, k]):
                sample_number = np.sum(M_array[0, :k]) + i
                member_ilo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel0_ilevel{k-1}_sample{sample_number}.nc', engine='netcdf4')
                member_ihi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel0_ilevel{k}_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[0, k-1]:
                    t_interpolated_member_ilo = interpolate.interpolate_time(reference_data['time'].values, member_ilo, output_vars=output_vars)
                else:
                    t_interpolated_member_ilo = member_ilo
                if space_interpolation[0, k-1]:
                    interpolated_member_ilo = interpolate.interpolate_space(reference_data, t_interpolated_member_ilo, output_vars=output_vars)
                else:
                    interpolated_member_ilo = t_interpolated_member_ilo
                
                if time_interpolation[0, k]:
                    t_interpolated_member_ihi = interpolate.interpolate_time(reference_data['time'].values, member_ihi, output_vars=output_vars)
                else:
                    t_interpolated_member_ihi = member_ihi
                if space_interpolation[0, k]:
                    interpolated_member_ihi = interpolate.interpolate_space(reference_data, t_interpolated_member_ihi, output_vars=output_vars)
                else:
                    interpolated_member_ihi = t_interpolated_member_ihi

                # add to ensemble mean

                for variable in output_vars:
                    ensemble_mean[variable].values += (interpolated_member_ihi[variable].values - interpolated_member_ilo[variable].values) / M_array[0, k]
                    if squared:
                        ensemble_mean_of_square[variable].values += (np.power(interpolated_member_ihi[variable].values, 2) - \
                                                                np.power(interpolated_member_ilo[variable].values, 2)) / M_array[0, k]

        ## HIGHER OUTER LEVEL CORRECTIONS ##

        for l in range(1, num_outer_levels):
            # Inner level 0
            for i in range(M_array[l, 0]):
                sample_number = np.sum(M_array[:l, :]) + i
                member_olo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l-1}_ilevel0_sample{sample_number}.nc', engine='netcdf4')
                member_ohi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l}_ilevel0_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[l-1, 0]:
                    t_interpolated_member_olo = interpolate.interpolate_time(reference_data['time'].values, member_olo, output_vars=output_vars)
                else:
                    t_interpolated_member_olo = member_olo
                if space_interpolation[l-1, 0]:
                    interpolated_member_olo = interpolate.interpolate_space(reference_data, t_interpolated_member_olo, output_vars=output_vars)
                else:
                    interpolated_member_olo = t_interpolated_member_olo

                if time_interpolation[l, 0]:
                    t_interpolated_member_ohi = interpolate.interpolate_time(reference_data['time'].values, member_ohi, output_vars=output_vars)
                else:
                    t_interpolated_member_olo = member_olo
                if space_interpolation[l, 0]:
                    interpolated_member_ohi = interpolate.interpolate_space(reference_data, t_interpolated_member_ohi, output_vars=output_vars)
                else:
                    interpolated_member_ohi = t_interpolated_member_ohi

                # add to ensemble mean

                for variable in output_vars:
                    ensemble_mean[variable].values += (interpolated_member_ohi[variable].values - interpolated_member_olo[variable].values) / M_array[l, 0]
                    if squared:
                        ensemble_mean_of_square[variable].values += (np.power(interpolated_member_ohi[variable].values, 2) - \
                                                                np.power(interpolated_member_olo[variable].values, 2)) / M_array[l, 0]
            
            # Higher inner level corrections
            for k in range(1, num_inner_levels):
                for i in range(M_array[l, k]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :k]) + i

                    member_olo_ilo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l-1}_ilevel{k-1}_sample{sample_number}.nc', engine='netcdf4')
                    member_olo_ihi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l-1}_ilevel{k}_sample{sample_number}.nc', engine='netcdf4')
                    member_ohi_ilo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l}_ilevel{k-1}_sample{sample_number}.nc', engine='netcdf4')
                    member_ohi_ihi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l}_ilevel{k}_sample{sample_number}.nc', engine='netcdf4')

                    # interpolate if necessary
                    if time_interpolation[l-1, k-1]:
                        t_interpolated_member_olo_ilo = interpolate.interpolate_time(reference_data['time'].values, member_olo_ilo, output_vars=output_vars)
                    else:
                        t_interpolated_member_olo_ilo = member_olo_ilo
                    if space_interpolation[l-1, k-1]:
                        interpolated_member_olo_ilo = interpolate.interpolate_space(reference_data, t_interpolated_member_olo_ilo, output_vars=output_vars)
                    else:
                        interpolated_member_olo_ilo = t_interpolated_member_olo_ilo
                    
                    if time_interpolation[l-1, k]:
                        t_interpolated_member_olo_ihi = interpolate.interpolate_time(reference_data['time'].values, member_olo_ihi, output_vars=output_vars)
                    else:
                        t_interpolated_member_olo_ihi = member_olo_ihi
                    if space_interpolation[l-1, k]:
                        interpolated_member_olo_ihi = interpolate.interpolate_space(reference_data, t_interpolated_member_olo_ihi, output_vars=output_vars)
                    else:
                        interpolated_member_olo_ihi = t_interpolated_member_olo_ihi

                    if time_interpolation[l, k-1]:
                        t_interpolated_member_ohi_ilo = interpolate.interpolate_time(reference_data['time'].values, member_ohi_ilo, output_vars=output_vars)
                    else:
                        t_interpolated_member_ohi_ilo = member_ohi_ilo
                    if space_interpolation[l, k-1]:
                        interpolated_member_ohi_ilo = interpolate.interpolate_space(reference_data, t_interpolated_member_ohi_ilo, output_vars=output_vars)
                    else:
                        interpolated_member_ohi_ilo = t_interpolated_member_ohi_ilo

                    if time_interpolation[l, k]:
                        t_interpolated_member_ohi_ihi = interpolate.interpolate_time(reference_data['time'].values, member_ohi_ihi, output_vars=output_vars)
                    else:
                        t_interpolated_member_ohi_ihi = member_ohi_ihi
                    if space_interpolation[l, k]:
                        interpolated_member_ohi_ihi = interpolate.interpolate_space(reference_data, t_interpolated_member_ohi_ihi, output_vars=output_vars)
                    else:
                        interpolated_member_ohi_ihi = t_interpolated_member_ohi_ihi

                    # add to ensemble mean

                    for variable in output_vars:
                        ensemble_mean[variable].values += (interpolated_member_ohi_ihi[variable].values - interpolated_member_olo_ihi[variable].values - \
                                                      interpolated_member_ohi_ilo[variable].values + interpolated_member_olo_ilo[variable].values) / M_array[l, k]
                        if squared:
                            ensemble_mean_of_square[variable].values += (np.power(interpolated_member_ohi_ihi[variable].values, 2) - \
                                                                    np.power(interpolated_member_olo_ihi[variable].values, 2) - \
                                                                    np.power(interpolated_member_ohi_ilo[variable].values, 2) + \
                                                                    np.power(interpolated_member_olo_ilo[variable].values, 2)) / M_array[l, k]
        if squared:
            return ensemble_mean, ensemble_mean_of_square
        else:
            return ensemble_mean

    
    @staticmethod
    def NML_representative_ensemble(path, run_name, M_array, reference_data, time_interpolation, space_interpolation, longitude_point, latitude_point, time_index, output_var, size):
        """Generates approximate high-level (full model) independent and identically distributed ensemble at a single point in time and space,
         following a method analogous to the method by Gregory and Cotter (2017).
          
        - path (str):                       path to folder containing all the EPS runs;
        - run_name(str):                    name of the EPS run of interest;
        - reference_data (xArray dataset):  dataset containing fine grids to facilitate interpolation;
        - time_interpolation (np.ndarray):  array of booleans indicating whether approximations should be interpolated in time. Row index indicates 
                                            outer approximation level and column index indicates inner level;
        - space_interpolation (np.ndarray): array of booleans indicating whether approximations should be interpolated in space. Row index
                                            indicates outer approximation level and column index indicates inner level;
        - longitude_point (float):          longitude of the point of interest;
        - latitude_point (float):           latitude of the point of interest;
        - time_index (int):                 time index of interest;
        - output_var (str):                 output variable of interest;
        - size (int):                       size of the representative ensemble. 
        """   

        # Compute latitude and longitude indices
        lat_index = np.argmin(np.absolute(np.full(reference_data['latitude'].values.shape, latitude_point) - reference_data['latitude'].values))
        lon_index = np.argmin(np.absolute(np.full(reference_data['longitude'].values.shape, longitude_point) - reference_data['longitude'].values))
        point_of_interest = (time_index, lat_index, lon_index)

        num_outer_levels = M_array.shape[0]
        num_inner_levels = M_array.shape[1]

        # Gather data
        point_ensemble = np.zeros(M_array[0,0] + 2 * np.sum(M_array[0,1:]) + 2 * np.sum(M_array[1:,0]) + 4 * np.sum(M_array[1:,1:]))

        ## OUTER LEVEL 0 ##
        for i in range(M_array[0, 0]):
            sample_number = i
            point_ensemble_index = i
            member = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel0_ilevel0_sample{sample_number}.nc', engine='netcdf4')

            # interpolate if necessary
            if time_interpolation[0, 0]:
                t_interpolated_member = interpolate.interpolate_time(reference_data['time'].values, member, output_vars=[output_var])
            else:
                t_interpolated_member = member
            if space_interpolation[0, 0]:
                interpolated_member = interpolate.interpolate_space(reference_data, t_interpolated_member, output_vars=[output_var])
            else:
                interpolated_member = t_interpolated_member
            
            # add to point ensemble
            point_ensemble[point_ensemble_index] = interpolated_member[output_var].values[point_of_interest]

        for k in range(1, num_inner_levels):
            for i in range(M_array[0, k]):
                sample_number = np.sum(M_array[0, :k]) + i
                point_ensemble_index_ilo = M_array[0,0] + 2 * np.sum(M_array[0, :k]) + i
                point_ensemble_index_ihi = point_ensemble_index_ilo + M_array[0, k] + i

                member_ilo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel0_ilevel{k-1}_sample{sample_number}.nc', engine='netcdf4')
                member_ihi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel0_ilevel{k}_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary

                if time_interpolation[0, k-1]:
                    t_interpolated_member_ilo = interpolate.interpolate_time(reference_data['time'].values, member_ilo, output_vars=[output_var])
                else:
                    t_interpolated_member_ilo = member_ilo
                if space_interpolation[0, k-1]:
                    interpolated_member_ilo = interpolate.interpolate_space(reference_data, t_interpolated_member_ilo, output_vars=[output_var])
                else:
                    interpolated_member_ilo = t_interpolated_member_ilo
                
                if time_interpolation[0, k]:
                    t_interpolated_member_ihi = interpolate.interpolate_time(reference_data['time'].values, member_ihi, output_vars=[output_var])
                else:
                    t_interpolated_member_ihi = member_ihi
                if space_interpolation[0, k]:
                    interpolated_member_ihi = interpolate.interpolate_space(reference_data, t_interpolated_member_ihi, output_vars=[output_var])
                else:
                    interpolated_member_ihi = t_interpolated_member_ihi

                # add to point ensemble

                point_ensemble[point_ensemble_index_ilo] = interpolated_member_ilo[output_var].values[point_of_interest]
                point_ensemble[point_ensemble_index_ihi] = interpolated_member_ihi[output_var].values[point_of_interest]

        ## HIGHER OUTER LEVEL CORRECTIONS ##
        for l in range(1, num_outer_levels):
            for i in range(M_array[l, 0]):
                sample_number = np.sum(M_array[:l, :]) + i
                point_ensemble_index_olo = M_array[0,0] + 2 * np.sum(M_array[0, 1:]) + 2 * np.sum(M_array[1:l, 0]) + 4 * np.sum(M_array[1:l, 1:]) + i
                point_ensemble_index_ohi = point_ensemble_index_olo + M_array[l, 0]

                member_olo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l-1}_ilevel0_sample{sample_number}.nc', engine='netcdf4')
                member_ohi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l}_ilevel0_sample{sample_number}.nc', engine='netcdf4')

                # interpolate if necessary
                if time_interpolation[l-1, 0]:
                    t_interpolated_member_olo = interpolate.interpolate_time(reference_data['time'].values, member_olo, output_vars=[output_var])
                else:
                    t_interpolated_member_olo = member_olo
                if space_interpolation[l-1, 0]:
                    interpolated_member_olo = interpolate.interpolate_space(reference_data, t_interpolated_member_olo, output_vars=[output_var])
                else:
                    interpolated_member_olo = t_interpolated_member_olo

                if time_interpolation[l, 0]:
                    t_interpolated_member_ohi = interpolate.interpolate_time(reference_data['time'].values, member_ohi, output_vars=[output_var])
                else:
                    t_interpolated_member_ohi = member_ohi
                if space_interpolation[l, 0]:
                    interpolated_member_ohi = interpolate.interpolate_space(reference_data, t_interpolated_member_ohi, output_vars=[output_var])
                else:
                    interpolated_member_ohi = t_interpolated_member_ohi

                # add to point ensemble
                point_ensemble[point_ensemble_index_olo] = interpolated_member_olo[output_var].values[point_of_interest]
                point_ensemble[point_ensemble_index_ohi] = interpolated_member_ohi[output_var].values[point_of_interest]

            for k in range(1, num_inner_levels):
                for i in range(M_array[l, k]):
                    sample_number = np.sum(M_array[:l, :]) + np.sum(M_array[l, :k]) + i
                    start_index = M_array[0,0] + 2 * np.sum(M_array[0, 1:k]) + 2 * np.sum(M_array[1:l, 0]) + 4 * np.sum(M_array[1:l, 1:k]) # here is where this outer level starts
                    point_ensemble_index_olo_ilo = start_index + 2 * M_array[l, 0] + 4*np.sum(M_array[l, 1:k]) + i
                    point_ensemble_index_ohi_ilo = point_ensemble_index_olo_ilo + M_array[l, k]
                    point_ensemble_index_olo_ihi = point_ensemble_index_ohi_ilo + M_array[l, k]
                    point_ensemble_index_ohi_ihi = point_ensemble_index_olo_ihi + M_array[l, k]

                    member_olo_ilo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l-1}_ilevel{k-1}_sample{sample_number}.nc', engine='netcdf4')
                    member_olo_ihi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l-1}_ilevel{k}_sample{sample_number}.nc', engine='netcdf4')
                    member_ohi_ilo = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l}_ilevel{k-1}_sample{sample_number}.nc', engine='netcdf4')
                    member_ohi_ihi = xr.open_dataset(path + f'{run_name}\\output\\{run_name}_olevel{l}_ilevel{k}_sample{sample_number}.nc', engine='netcdf4')

                    # interpolate if necessary

                    if time_interpolation[l-1, k-1]:
                        t_interpolated_member_olo_ilo = interpolate.interpolate_time(reference_data['time'].values, member_olo_ilo, output_vars=[output_var])
                    else:
                        t_interpolated_member_olo_ilo = member_olo_ilo
                    if space_interpolation[l-1, k-1]:
                        interpolated_member_olo_ilo = interpolate.interpolate_space(reference_data, t_interpolated_member_olo_ilo, output_vars=[output_var])
                    else:
                        interpolated_member_olo_ilo = t_interpolated_member_olo_ilo

                    if time_interpolation[l-1, k]:
                        t_interpolated_member_olo_ihi = interpolate.interpolate_time(reference_data['time'].values, member_olo_ihi, output_vars=[output_var])
                    else:
                        t_interpolated_member_olo_ihi = member_olo_ihi
                    if space_interpolation[l-1, k]:
                        interpolated_member_olo_ihi = interpolate.interpolate_space(reference_data, t_interpolated_member_olo_ihi, output_vars=[output_var])
                    else:
                        interpolated_member_olo_ihi = t_interpolated_member_olo_ihi

                    if time_interpolation[l, k-1]:
                        t_interpolated_member_ohi_ilo = interpolate.interpolate_time(reference_data['time'].values, member_ohi_ilo, output_vars=[output_var])
                    else:
                        t_interpolated_member_ohi_ilo = member_ohi_ilo
                    if space_interpolation[l, k-1]:
                        interpolated_member_ohi_ilo = interpolate.interpolate_space(reference_data, t_interpolated_member_ohi_ilo, output_vars=[output_var])
                    else:
                        interpolated_member_ohi_ilo = t_interpolated_member_ohi_ilo

                    if time_interpolation[l, k]:
                        t_interpolated_member_ohi_ihi = interpolate.interpolate_time(reference_data['time'].values, member_ohi_ihi, output_vars=[output_var])
                    else:
                        t_interpolated_member_ohi_ihi = member_ohi_ihi
                    if space_interpolation[l, k]:
                        interpolated_member_ohi_ihi = interpolate.interpolate_space(reference_data, t_interpolated_member_ohi_ihi, output_vars=[output_var])
                    else:
                        interpolated_member_ohi_ihi = t_interpolated_member_ohi_ihi

                    # add to point ensemble
                    point_ensemble[point_ensemble_index_olo_ilo] = interpolated_member_olo_ilo[output_var].values[point_of_interest]
                    point_ensemble[point_ensemble_index_ohi_ilo] = interpolated_member_ohi_ilo[output_var].values[point_of_interest]
                    point_ensemble[point_ensemble_index_olo_ihi] = interpolated_member_olo_ihi[output_var].values[point_of_interest]
                    point_ensemble[point_ensemble_index_ohi_ihi] = interpolated_member_ohi_ihi[output_var].values[point_of_interest]

        # Sort the point ensemble in the correct way
        sorted_ensemble = np.zeros(point_ensemble.shape)

        sorted_ensemble[:M_array[0,0]] = np.sort(point_ensemble[:M_array[0,0]])

        for k in range(1, num_inner_levels):
            start_index = M_array[0,0] + 2 * np.sum(M_array[0, 1:k])
            sorted_ensemble[start_index:(start_index + M_array[0, k])] = np.sort(point_ensemble[start_index:(start_index + M_array[0,k])])
            sorted_ensemble[(start_index + M_array[0,k]):(start_index + 2*M_array[0,k])] = np.sort(point_ensemble[(start_index+M_array[0,k]):(start_index + 2*M_array[0,k])]) 

        for l in range(1, num_outer_levels):
            start_index = M_array[0,0] + 2 * np.sum(M_array[1:l, 0]) + 2 * np.sum(M_array[0, 1:]) + 4 * np.sum(M_array[1:l, 1:])
            sorted_ensemble[start_index:(start_index + M_array[l,0])] = np.sort(point_ensemble[start_index:(start_index + M_array[l,0])])
            sorted_ensemble[(start_index + M_array[l,0]):(start_index+2*M_array[l,0])] = np.sort(point_ensemble[(start_index + M_array[l,0]):(start_index+2*M_array[l,0])])

            for k in range(1, num_inner_levels):
                start_index = M_array[0,0] + 2 * np.sum(M_array[1:l, 0]) + 2 * np.sum(M_array[0, 1:]) + 4 * np.sum(M_array[1:l, 1:]) + \
                    2 * M_array[l, 0] + 4 * M_array[l, 1:k]
                sorted_ensemble[start_index:(start_index + M_array[l,k])] = np.sort(point_ensemble[start_index:(start_index + M_array[l,k])])
                sorted_ensemble[(start_index + M_array[l,k]):(start_index + 2*M_array[l,k])] = np.sort(point_ensemble[(start_index + M_array[l,k]):(start_index + 2*M_array[l,k])])
                sorted_ensemble[(start_index + 2*M_array[l,k]):(start_index + 3*M_array[l,k])] = np.sort(point_ensemble[(start_index + 2*M_array[l,k]):(start_index + 3*M_array[l,k])])
                sorted_ensemble[(start_index + 3*M_array[l,k]):(start_index + 4*M_array[l,k])] = np.sort(point_ensemble[(start_index + 3*M_array[l,k]):(start_index + 4*M_array[l,k])])
                
                                                                                                              
        # Generate representative ensemble

        repr_ensemble = np.zeros(size)

        for n in range(size):
            u = np.random.uniform()

        repr_ensemble[n] += sorted_ensemble[int(np.floor(M_array[0,0] * u))]

        for k in range(1, num_inner_levels):
            start_index = M_array[0,0] + 2 * np.sum(M_array[0, 1:k])
            repr_ensemble[n] += -sorted_ensemble[start_index + int(np.floor(M_array[0,k] * u))] + \
                                sorted_ensemble[start_index + M_array[0,k] + int(np.floor(M_array[0,k] * u))]
            
        for l in range(1, num_outer_levels):
            start_index = M_array[0,0] + 2 * np.sum(M_array[1:l, 0]) + 2 * np.sum(M_array[0, 1:]) + 4 * np.sum(M_array[1:l, 1:])
            repr_ensemble[n] += -sorted_ensemble[start_index + int(np.floor(M_array[l,0] * u))] + \
                                sorted_ensemble[start_index + M_array[l, 0] + int(np.floor(M_array[l,0] * u))]

            for k in range(1, num_inner_levels):
                start_index = M_array[0,0] + 2 * np.sum(M_array[1:l, 0]) + 2 * np.sum(M_array[0, 1:]) + 4 * np.sum(M_array[1:l, 1:]) + \
                    2 * M_array[l, 0] + 4 * M_array[l, 1:k]
                
                repr_ensemble[n] += sorted_ensemble[start_index + 3*M_array[l,k] + int(np.floor(M_array[l,k] * u))] - \
                                    sorted_ensemble[start_index + 2*M_array[l,k] + int(np.floor(M_array[l,k] * u))] - \
                                    sorted_ensemble[start_index + M_array[l,k] + int(np.floor(M_array[l,k] * u))] + \
                                    sorted_ensemble[start_index + int(np.floor(M_array[l,k] * u))]
                
        return repr_ensemble

            
        
##############
## OLD CODE ##
##############

# def ensemble_mean_map(ensemble, timestep, parameter):
#     """Computes the ensemble mean of a parameter from an independent, identically distributed ensemble at a single timestep.
#     Does not work for multilevel ensembles.
    
#     - ensemble (list):      list of xArray datasets forming an independent, identically distributed ensemble of the output distribution;
#     - timestep (int):       time index the ensemble members should be evaluated on;
#     - parameter (str):      name of the output parameter of interest.
#     """

#     summed_ensemble = np.zeros(ensemble[0][parameter][timestep, :, :].values.shape)

#     for i in range(len(ensemble)):
#         summed_ensemble += ensemble[i][parameter][timestep, :, :].values

#     return summed_ensemble / len(ensemble)


# def ensemble_spread_map(ensemble, timestep, parameter, mean=None):
#     """Computes the ensemble spread of a parameter from an independent, identically distributed ensemble at a single timestep.
#     Does not work for multilevel ensembles."""

#     if mean is None:
#         ensemble_mean = ensemble_mean_map(ensemble, timestep, parameter)
#     else:
#         ensemble_mean = mean[:]

#     square_summed_ensemble = np.zeros(ensemble_mean.shape)

#     for i in range(len(ensemble)):
#         square_summed_ensemble += np.power(ensemble[i][parameter][timestep, :, :].values - ensemble_mean, 2)

#     return square_summed_ensemble / len(ensemble)


########################
## DYSFUNCTIONAL CODE ##
########################

##################################################################################################################################
# ## SHOW ANIMATION MAP (DOESN'T WORK)##

# def animate_map(lon, lat, time, map, cmap='viridis', title=None, label=None, quiver=None, stride=None, savestr=None):
#     fig, ax = plt.subplots(figsize = (10,6), subplot_kw={'projection': ccrs.PlateCarree()})

#     global_max = np.amax(map)
#     global_min = np.amin(map)

#     ax.gridlines(draw_labels=True)
#     ax.coastlines()

#     plot_colormap(fig, ax, lon, lat, map[0,:,:], title=title, quantitystr=label, vmin=global_min, vmax=global_max, cmap=cmap)
#     if (not quiver is None) and (not stride is None):
#         plot_arrows_on_colormap(ax, lon, lat, stride, quiver)

#     def frame(i):
#         ax.clear()
#         plot_colormap(fig, ax, lon, lat, map[i,:,:], title=title, quantitystr=label, vmin=global_min, vmax=global_max, cmap=cmap)
#         ax.gridlines(draw_labels=True)
#         ax.coastlines()
#         if quiver and stride:
#             plot_arrows_on_colormap(ax, lon, lat, stride, quiver)

#     ani = anim.FuncAnimation(fig, frame, frames=time.values.shape[0])
#     plt.show()
##########################################################################################################################


    

    

