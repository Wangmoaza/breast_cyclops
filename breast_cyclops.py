import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

############ TODO ##############
# 95% lowerbound upperbound
#
#
################################

def find_peak(values):
    """ Returns the indices of the boundaries of peaks. A peak is defined as
    local maxima where the height of the local maxima to its surrounding 
    local minima are greater than (global max - global min)/10.

    args:
        values (numpy 1D array) : count values
    
    returns:
        2D array where each element contains [lower boundary, upper boundary] of the peak.
        The boundaries are indices not the value themselves.
    """

    # values should not contain NaN values
    x = values
    # argrelextrema returns indices of local maxima/minima
    # local maxima
    local_max = argrelextrema(x, np.greater)
    local_max = local_max[0]
    # local minima
    local_min = argrelextrema(x, np.less)
    local_min = local_min[0]

    lower = 0
    upper = 0
    peaks = []
    # for each local maxima, if height of the local maxima to its surrounding local minima is
    # greater than (global max - global min)/10, consider it as a peak
    for maximum in local_max:
        while local_min[lower] < maximum:
            lower += 1
        upper = lower
        lower -= 1
        while (upper < local_min.shape[0]) and (local_min[upper] < maximum):
            upper += 1
        # FIXME height of peak > (global max - global min)/10
        if x[maximum] - min(x[local_min[lower]], x[local_min[upper]]) > (x.max() - x.min()) / 10: 
            peaks.append(np.array([local_min[lower], local_min[upper]]))
    ### END - for maximum
    return np.array(peaks)
### END - find_peak


def main():
    cnv_df = pd.read_table('../CCLE_BR_lines_CNV.txt', sep='\t', header=0, index_col=0)
    hist, bin_edges = np.histogram(cnv['BT20_BREAST'].values, bins=100)
    sample_hist = pd.DataFrame(data=hist, index=bin_edges[:-1], columns=['count'])
    # 5-bin moving average
    smooth_hist = sample_hist.rolling(window=5).mean()
    smooth_hist = smooth_hist[smooth_hist['count'].notnull()] # remove NaNs
    smooth_hist.reset_index(inplace=True)
    smooth_hst.columns = ['lower_boound', 'count']
    counts = smooth_hist['count'].values
    find_peak(counts)


