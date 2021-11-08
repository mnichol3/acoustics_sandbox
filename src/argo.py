from argopy import DataFetcher as ArgoDataFetcher
import numpy as np


def get_profile(bbox, dtg_start, dtg_end, min_db=0.0, max_db=2000.0):
    '''
    Get a vertical profile of temperature, salinity, and pressure from an
    Argo float.

    Parameters
    -----------
    bbox : list of float
        List of floats defining a geographical bounding box to select Argo float
        samples from. Format: [min_lon, max_lon, min_lat, max_lat].
    dtg_start : str
        8-digit datetime group definining the start of the Argo sample period.
        Format: YYYY-MM-DD.
    dtg_end : str
        8-digit datetime group defining the end of the Argo sample period.
        Format: YYYY-MM-DD.
    min_db : float, optional
        Minimum pressure to include in Argo samples. Default is 0.
    max_db : float, optional
        Maximum float to include in Argo samples. Default is 1500.

    Returns
    --------
    xarray.DataArray
    '''
    argo_params = bbox.append(min_db, max_db, dtg_start, dtg_end)

    ## Get the argo data for our specified space & time
    argo_loader = ArgoDataFetcher().region(argo_params)
    argo_pts = argo_loader.to_xarray()
    argo_profiles = argo_pts.argo.point2profile()

    ## Select a profile
    for argo_profiler in list(argo_profiles.N_PROF.data):
        curr_profile = argo_profiles.isel(N_PROF=argo_profiler)
        ## Make sure we have a vertical profile with ample data points
        curr_press = curr_profile.PRES.data
        num_samples = curr_pres[~np.isnan(curr_pres)]
        if num_samples >= max_db / 2:
            break

    return argo_profile
