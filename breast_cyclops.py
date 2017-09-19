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
        2D array where each element contains [lower bound, summit, upper bound] of the peak.
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
    # FIXME height threshold of peak (global max - global min)/10
    height_threshold = (x.max() - x.min())/10
    # for each local maxima, if height of the local maxima to its surrounding local minima is
    # greater than (global max - global min)/10, consider it as a peak
    for maximum in local_max:
        while local_min[lower] < maximum:
            lower += 1
        upper = lower
        lower -= 1
        while (upper < local_min.shape[0]) and (local_min[upper] < maximum):
            upper += 1
        if x[maximum] - min(x[local_min[lower]], x[local_min[upper]]) > height_threshold: 
            peaks.append(np.array([local_min[lower], maximum, local_min[upper]]))
    ### END - for maximum
    return np.array(peaks)
### END - find_peak


def _get95bound(smooth_hist, low_idx, high_idx, direction='upper'):
    # FIXME What does 95% bounds mean? Possibly remove this
    # 5% of counts in the peak range
    target = smooth_hist.loc[low_idx:high_idx, 'count'].sum() * 0.05
    subtotal = 0
    cutoff = 0
    a, b, d = 0, 0, 0
    if direction == 'lower':
        a, b, d = low_idx, high_idx + 1, 1

    elif direction == 'upper':
        a, b, d = high_idx, low_idx - 1, -1

    else:
        print "Error: Invalid direction."
        return

    for i in range(a, b, d):
        subtotal += smooth_hist.loc[i, 'count']
        if subtotal > target:
            x1, x2 = smooth_hist.loc[i-1, 'lower_bound'], smooth_hist.loc[i, 'lower_bound']
            y1, y2 = smooth_hist.loc[i-1, 'count'], smooth_hist[i, 'count']
            cutoff = (x2 - x1)/(y2 - y1) * (target - y1) + x1
            return cutoff
        ### END - if
    ### END - for i
### END - get95bound

def get95bound(smooth_hist, low_idx, high_idx, direction='upper'):
    # intergrate the counts in [low_idx, high_idx] using trapeziodal rule
    target = np.trapz(smooth_hist[low_idx:high_idx, 'count'].values, 
                      x=smooth_hist[low_idx:high_idx, 'lower_bound'].values) * 0.05



def main():
    cnv_df = pd.read_table('../CCLE_BR_lines_CNV.txt', sep='\t', header=0, index_col=0)
    hist, bin_edges = np.histogram(cnv['BT20_BREAST'].values, bins=100)
    sample_hist = pd.DataFrame(data=hist, index=bin_edges[:-1], columns=['count'])
    # 5-bin moving average
    smooth_hist = sample_hist.rolling(window=5).mean()
    smooth_hist = smooth_hist[smooth_hist['count'].notnull()] # remove NaNs
    smooth_hist.reset_index(inplace=True)
    smooth_hist.columns = ['lower_bound', 'count']
    counts = smooth_hist['count'].values
    
    peak_indices = find_peak(counts)
    peak_loss, peak_neutral = np.nan, np.nan
    for peak in peak_indices:
        value = smoothe_hist.loc[peak[1], 'lower_bound'] 
        if -0.4 < value and value < -0.05:
            peak_loss = peak
        if -0.05 < value and value < 0.05:
            peak_neutral = peak

    if np.isfinite(peak_loss) and np.isfinite(peak_neutral):
        #TODO

    else:
        

