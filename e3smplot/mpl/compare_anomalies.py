#!/usr/bin/env python3
from e3smplot.e3sm_utils import open_dataset, get_common_time_range

def plot_anomalies(means, labels, **kwargs):
    '''
    Plot anomalies. Accept a list of means and labels. Make a scatter of the
    distribution of each mean
    '''
    return None

def plot_scatter(data, labels, **kwargs):
    return None
    
def field_in_dataset():
    return None

def main(inputfiles, plotname, *varnames, labels=None, **kwargs):

    # Open datasets; note that inputfiles may be a list of lists. In this case,
    # each item in inputfiles is a list of files for a specific case;
    # open_dataset should use xarray.open_mfdatasets() to concatenate these to a
    # single dataset
    datasets = [open_dataset(f) for f in inputfiles] 

    # Subset for overlapping time period
    t1, t2 = get_common_time_range(datasets)

    # Compute daily means

    # Compute global means

    # Compute anomalies

    # Make plots

if __name__ == '__main__':
    import plac; plac.call(main)
