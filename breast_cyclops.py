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

    lower, upper = 0, 0
    peaks = []
    lower_minima_idx, upper_minima_idx = -1, -1
    # FIXME height threshold of peak (global max - global min)/10
    height_threshold = (x.max() - x.min())/10

    # for each local maxima, if height of the local maxima to its surrounding local minima is
    # greater than (global max - global min)/10, consider it as a peak
    for maximum in local_max:
        try:
            while local_min[lower] < maximum:
                lower += 1
            upper = lower
            lower -= 1
            while (upper < local_min.shape[0]) and (local_min[upper] < maximum):
                upper += 1 
            lower_minima_idx = local_min[lower]
            upper_minima_idx = local_min[upper]

        except IndexError:
            print "index error"
            lower -= 1
            lower_minima_idx = local_min[lower]
            upper_minima_idx = x.shape[0] - 1
        
        if lower_minima_idx < 0 or upper_minima_idx < 0:
            print "Error: indices should be > 0"
            return
        # FIXME height of peak > (global max - global min)/10            
        if x[maximum] - min(x[lower_minima_idx], x[upper_minima_idx]) > height_threshold: 
            peaks.append(np.array([lower_minima_idx, maximum, upper_minima_idx]))

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
### END - _get95bound

def get95bound(smooth_hist, low_idx, high_idx, direction='upper'):
    """ Calculate upper/lower 95% bounds for given distribution.

        args:
            smooth_hist (pandas Series) : peak. indices are log2CN lowerbounds and values are counts
            low_idx (int) : index of low bound for the peak in smooth_hist
            high_idx (int) : index of high bound for the peak in smooth_hist
            direction (string) : (default='upper') upper or lower 95% bound

        returns:
            log2CN cutoff for upper/lower 95% bound
    """

    # intergrate the counts in [low_idx, high_idx] using trapeziodal rule
    target = np.trapz(smooth_hist.loc[low_idx:high_idx, 'count'].values, 
                      x=smooth_hist.loc[low_idx:high_idx, 'lower_bound'].values) * 0.05
    area, cutoff = 0, 0
    a, b, d = 0, 0, 0
    if direction == 'lower':
        a, b, d = low_idx, high_idx + 1, 1

    elif direction == 'upper':
        a, b, d = high_idx, low_idx - 1, -1

    else:
        print "Error: Invalid direction."
        return

    # FIXME Currently linear serach. Should change to binary search?
    # FIXME Need to check if this is correct...
    for i in range(a, b, d):
        if d > 0:
            new_area = np.trapz(smooth_hist.loc[i:i+d, 'count'], smooth_hist.loc[i:i+d, 'lower_bound'])
        else:
            new_area = np.trapz(smooth_hist.loc[i+d:i, 'count'], smooth_hist.loc[i+d:i, 'lower_bound'])
        target -= new_area
        if target < 0:
            target += new_area
            x1, x2 = smooth_hist.loc[i, 'lower_bound'], smooth_hist.loc[i+d, 'lower_bound']
            y1, y2 = smooth_hist.loc[i, 'count'], smooth_hist.loc[i+d, 'count']
            slope = float(y2-y1)/(x2-x1)
            discriminant = np.sqrt(y1**2 + (2 * d * slope * target))
            cutoff = (-y1 + (slope * x1) - discriminant) / slope
            if cutoff < x1:
                cutoff += 2 * (discriminant / slope)
            return cutoff
        ### END - if
        elif target == 0:
            return smooth_hist.loc[i+d, 'lower_bound']
    ### END - for i
### END - get95bound

# FIXME needs correction check
def classify_loss(cnv_df):
    # partial loss (1), homozygous loss (2), copy neutral (0), uncalled(NaN) info stored in loss_df
    loss_df = pd.DataFrame(data=np.nan, columns = cnv_df.columns, index=cnv_df.index)
    cutoff_homo = -1.28

    for cell_line in cnv_df.columns:
        print cell_line, '...'
        hist, bin_edges = np.histogram(cnv_df[cell_line].values, bins=100)       
        sample_hist = pd.DataFrame(data=hist, index=bin_edges[:-1], columns=['count'])
        smooth_hist = sample_hist.rolling(window=5).mean() # 5-bin moving average
        smooth_hist = smooth_hist[smooth_hist['count'].notnull()] # remove NaNs
        smooth_hist.reset_index(inplace=True)
        smooth_hist.columns = ['lower_bound', 'count']
        counts = smooth_hist['count'].values
        
        peak_indices = find_peak(counts)
        peak_loss, peak_neutral = None, None
        cutoff_loss, cutoff_neutral = np.nan, np.nan

        # check if peaks are in desired ranges
        # FIXME If there are multiple peaks in the region, should select highest peak?
        # currently choose rightmost peak for peak_loss, leftmost peak for peak_neutral
        for peak in peak_indices:
            summit = smooth_hist.loc[peak[1], 'lower_bound'] 
            if -0.4 < summit and summit < -0.05:
                peak_loss = peak
            if -0.05 < summit and summit < 0.05:
                peak_neutral = peak
                break
        ### END - for peak

        # log2CN < cutoff_loss are considered copy loss
        # log2CN > cutoff_neutral are considered copy neutral
        # loss peak and neutral peak are in desired ranges
        if peak_loss is not None and peak_neutral is not None:
            cutoff_loss    = get95bound(smooth_hist, peak_loss[0], peak_loss[2], direction='upper')
            cutoff_neutral = get95bound(smooth_hist, peak_neutral[0], peak_neutral[2], direction='lower')
        else:
            cutoff_loss    = -0.4
            cutoff_neutral = - 0.2
        ### END - if else
        
        # set loss_df according to cutoffs
        condition_loss    = (cutoff_homo < cnv_df[cell_line]) & (cnv_df[cell_line] < cutoff_loss)
        condition_neutral = cnv_df[cell_line] > cutoff_neutral
        condition_homo    = cnv_df[cell_line] <= cutoff_homo       
        loss_df.loc[cnv_df[condition_homo].index, cell_line]    = 2
        loss_df.loc[cnv_df[condition_loss].index, cell_line]    = 1
        loss_df.loc[cnv_df[condition_neutral].index, cell_line] = 0
    ### END - for cell_line
    return loss_df
### END - classify_loss

def main():
    cnv_df = pd.read_table('../CCLE_BR_lines_CNV.txt', sep='\t', header=0, index_col=0)
    cnv_df.columns = cnv_df.columns.str[:-7] # drop '_BREAST' part in cell_lines
    loss_df = classify_loss(cnv_df)
    loss_df.to_csv('CCLE_BR_lines_loss.tsv', sep='\t')

if __name__ == '__main__':
    main()
