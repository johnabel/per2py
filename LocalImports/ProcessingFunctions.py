'''
File containing functions useful for data analysis.
'''


from __future__ import division

# imports
import os
import numpy  as np
import scipy as sp
import pandas as pd
import spectrum
from matplotlib import gridspec
import matplotlib.pyplot as plt

from LocalImports import PlotOptions as plo
from LocalImports import Bioluminescence as blu
from LocalImports import DecayingSinusoid as dsin

def generate_filenames_dict(input_folder, input_file,
                        pull_from_imagej, input_ij_extension='.csv'):
    """
    Creates dict of input data
    """
    input_dict = {'input_folder':input_folder,
                  'data_type':input_file,
                  'pull_from_imagej':pull_from_imagej,
                  'input_ij_file':input_file,
                  'input_ij_extension':input_ij_extension}

    #inputs - autoset
    input_dict['input_data'] = input_folder+input_file+input_ij_extension
    data_type = input_file

    #outputs - general
    output_extension = '.csv'
    output_folder = input_folder+'analysis_output/'
    try:
        os.mkdir(output_folder)
    except OSError:
        pass

    # outputs - locations of raw data
    input_dict['raw_signal'] = output_folder +data_type+'_signal'+output_extension
    input_dict['raw_xy'] = output_folder +data_type+'_XY'+output_extension

    # outputs - processed data: hodrick-prescott detrend, and hp detrend plus eigensmoothing
    input_dict['output_detrend'] = output_folder+data_type+'_signal_detrend'+output_extension
    input_dict['output_detrend_smooth'] = output_folder+data_type+'_signal_detrend_denoise'+output_extension
    input_dict['output_detrend_smooth_xy'] = output_folder+data_type+'_XY'+'_detrend_denoise'+output_extension
    input_dict['output_pgram'] = output_folder+data_type+'_lombscargle'+output_extension
    input_dict['output_cosine'] = output_folder+data_type+'_cosine'+output_extension
    input_dict['output_phases'] = output_folder+data_type+'_cosine_phases'+output_extension
    input_dict['output_cosine_params'] = output_folder+data_type+'_oscillatory_params'+output_extension
    input_dict['output_zscore'] = output_folder+data_type+'_zscore'+output_extension
    return input_dict

def load_imagej_file(input_data, raw_signal, raw_xy, input_ij_extension=None):
    """
    Function to load imageJ files and save them in raw_signal and raw_xy
    locations. input_ij_extension can be set manually, defaults to last
    four characters of input_data.
    """
    print "Assembling data from imageJ recording... time: ",
    if input_ij_extension is None:
        input_ij_extension = input_data[-4:]
    timer = plo.laptimer()
    # load the file, it only has one sheet
    if input_ij_extension=='.csv':
        ij_sheet =  pd.read_csv(input_data)
    elif input_ij_extension=='.xls':
        ij_file = pd.read_excel(input_data, sheet_name=None)
        ij_sheet = ij_file[ij_file.keys()[0]]

    # x and y are self-evident, cell id is the track_id, mean_intensityXX
    # is the biolum we want
    useful_frames = ['TRACK_ID', 'POSITION_X', 'POSITION_Y', 'FRAME']

    # delete the rest
    orig_keys = np.copy(ij_sheet.keys())
    for kdx, key in enumerate(orig_keys):
        if key not in useful_frames:
            if key[:14]=='MEAN_INTENSITY':
                ij_sheet.rename(columns={key: 'MEAN_INTENSITY'}, inplace=True)
            else:
                del ij_sheet[orig_keys[kdx]]
    #


    # collect unique cell ID nos
    cell_ids = ij_sheet.TRACK_ID.unique()

    # set up array that has times and frames as left two cols for lum
    frames = np.sort(ij_sheet.FRAME.unique())
    biolum_data = np.nan*np.ones((len(frames), len(cell_ids)+2))
    biolum_data[:,0] = frames*0.25
    biolum_data[:,1] = frames

    # set up array that has times and frames as left three cols for xy
    xy_data = np.nan*np.ones((len(frames), len(cell_ids)*2+3))
    xy_data[:,0] = frames*0.25
    xy_data[:,1] = frames

    # assemble the data for each (time, frame, x, y, intensity)
    for cdx,cell_id in enumerate(cell_ids):
        # get the df where the cell is
        cell_df = ij_sheet[ij_sheet['TRACK_ID'].isin([cell_id])]
        # assemble!
        for fdx,frame in enumerate(frames):
            try:
                biolum_data[fdx,cdx+2] = \
                    cell_df[cell_df['FRAME'].isin([frame])].MEAN_INTENSITY.item()
                    #biolum
                xy_data[fdx,cdx*2+3] = \
                    cell_df[cell_df['FRAME'].isin([frame])].POSITION_X.item()#X
                xy_data[fdx,cdx*2+4] = \
                    cell_df[cell_df['FRAME'].isin([frame])].POSITION_Y.item()#Y
            except ValueError:
                # pass if there is no data there
                pass

    # create two output DFs
    biolum_df = pd.DataFrame(data=biolum_data,
        columns = ['TimesH', 'Frame']+list(cell_ids))
    xy_df = pd.DataFrame(data=xy_data,
      columns = ['TimesH', 'Frame', 'NoMissing']+list(np.sort(list(cell_ids)*2)))
    xy_df.columns = pd.MultiIndex.from_arrays([xy_df.columns,
        ['','','Place_Cursor']+['X','Y']*len(cell_ids)])

    # save them
    biolum_df.to_csv(raw_signal, index=False)
    xy_df.to_csv(raw_xy, index=False)

    print str(np.round(timer(),1)) +' s'

    del ij_sheet
    del biolum_df, biolum_data
    del xy_df, xy_data

def import_data(raw_signal, raw_xy):
    timer = plo.laptimer()

    print "Importing data... time: ",
    header = np.genfromtxt(raw_signal, delimiter=',')[0, :]
    raw_times = np.genfromtxt(raw_signal, delimiter=',')[1:, 0]
    raw_data = np.genfromtxt(raw_signal, delimiter=',')[1:, 2:]
    locations = np.genfromtxt(raw_xy, delimiter=',')[2:, 3:]

    # delete blank rows and columns
    raw_data = delete_blank_rows(raw_data)
    raw_data = delete_blank_columns(raw_data)
    raw_times = raw_times[:len(raw_data)]
    locations = delete_blank_rows(locations)
    locations = delete_blank_columns(locations)
    print str(np.round(timer(),1))+"s"
    return raw_times, raw_data, locations, header

def delete_blank_columns(dataset):
    """
    Deletes the columns consisting solely of nans in a dataset.
    """
    finished = False
    column_idx = 1
    total_cols = len(dataset[0,:])-1
    while finished is False and column_idx < total_cols:
        # count backwards to delete blank columns, each time look at the last one
        col = dataset[:,-1]
        is_nan, _ = blu.nan_helper(col)
        if all(is_nan):
            dataset = dataset[:,:-1]
            column_idx+=1
        else:
            finished = True
    return dataset

def delete_blank_rows(dataset):
    """
    Deletes the columns consisting solely of nans in a dataset.
    """
    finished = False
    row_idx = 1
    total_rows = len(dataset[:,0])-1
    while finished is False and row_idx < total_rows:
        # count backwards to delete blank columns, each time look at the last one
        row = dataset[-1,:]
        is_nan, _ = blu.nan_helper(row)
        if all(is_nan):
            dataset = dataset[:-1,:]
            row_idx+=1
        else:
            finished = True
    return dataset

def alignment(original, denoised, d=40, dstart=0):
    """
    The eigensmoothing as written truncates some front-back data, as the
    input data input data is extended. This function realigns the data

    The sign is sometimes flipped on the reconstructed signal.
    This function figures out +/- alignment as well.

    dstart (default=0) tells where to start checking for alignment. Some
    of the data has an artifact, so later datapoints are more reliable
    for finding an alignment (otherwise the artifact will dominate).
    """
    original = original[dstart:]
    denoised = denoised[dstart:]
    errs = np.zeros(d-1)
    for idx in range(d-1):
        errs[idx] = np.linalg.norm(original-denoised[idx:-(d-idx-1)])
    errs_neg = np.zeros(d-1)
    for idx in range(d-1):
        errs_neg[idx] = np.linalg.norm(original+denoised[idx:-(d-idx-1)])
    pos_min = np.min(errs)
    neg_min = np.min(errs_neg)
    if neg_min < pos_min:
        return np.argmin(errs_neg), -1
    else:
        return np.argmin(errs), 1

def truncate_and_interpolate(times, data, locations, truncate_t=12,
                                outlier_std=5):
    """
    Removes the first "truncate_t" h artifact and interpolates
    any missing values and rejects outliers.
    """
    # truncate the data
    start = np.argmax(times-times[0]>=truncate_t)
    times = times[start:]
    data = data[start:]
    locations = locations[start:]

    # first, find where nans are in the data
    h1, h2 = blu.nan_helper(data)

    # create output data
    outdata = np.nan*np.ones(data.shape)

    # interpolate the intermediate points
    for i in range(len(data[0,:])):
        try:
            init_goodval = np.min(np.where(~h1[:,i])[0])
            end_goodval  = np.max(np.where(~h1[:,i])[0])
            # only interpolate where we have good values
            y = np.copy(data[init_goodval:end_goodval, i])

            # find outliers - ignore where there's nans
            #  and ignore the warning here
            with np.errstate(invalid='ignore'):
                outliers = np.abs(y - np.nanmean(y)) > outlier_std * np.nanstd(y)
            outlier_idx = np.where(outliers)[0]
            # if any exist, replace them
            if len(outlier_idx)>0:
                # put in nans
                for oi in outlier_idx:
                    y[oi] = np.nan

            # find nans
            nans, x= blu.nan_helper(y)
            # interpolate
            y[nans]= np.interp(x(nans), x(~nans), y[~nans])
            # replace old values with new
            outdata[init_goodval:end_goodval, i] = y
        except ValueError:
            # if there are no nans
            pass

    # should we trim the data here as well? we can just leave it
    return times, outdata, locations

def hp_detrend(times, interpolated_data):
    """
    Does a Hodrick-Prescott detrending always using an estimated period
    of 24h.
    """
    timer = plo.laptimer()
    print "Detrending data... time: ",
    detrended_data = np.zeros(interpolated_data.shape)
    trendlines = np.zeros(interpolated_data.shape)
    for idx,idata in enumerate(interpolated_data.T):
        # copy data for editing
        cell = np.copy(idata)
        model = blu.Bioluminescence(times, cell, period_guess=24)
        # all we care about is the HP filter, no need to fit everything
        model.filter()
        model.detrend()
        baseline = model.yvals['mean']

        # find where recording starts and there are no nans
        valid = ~np.isnan(cell)
        # subtrect baseline from recording
        cell[valid] = cell[valid] - baseline
        # put the detrended cell in the detrended_data matrix,
        # and same for baseline
        detrended_data[:,idx] = cell
        trendlines[valid,idx] = baseline
    print str(np.round(timer(),1))+"s"
    return times, detrended_data, trendlines

def butterworth_lowpass(times, data, cutoff_period=4):
    """
    Does a Butterworth low pass filtering of the data.
    """
    timer = plo.laptimer()
    print "Butterworth filter... time: ",
    denoised_data = np.zeros(data.shape)
    for idx in range(len(data.T)):
        # copy data for editing
        idata = np.copy(data[:,idx])
        valid = ~np.isnan(idata)
        filtts, filtdata = blu.lowpass_filter(times[valid], idata[valid],
                                        cutoff_period=cutoff_period)

        # subtrect baseline from recording
        denoised_data[valid,idx] = filtdata

    print str(np.round(timer(),1))+"s"
    return times, denoised_data

def LS_pgram(times, ls_data, circ_low=18, circ_high=30, alpha=0.05):
    """Calculates a LS periodogram for each data sequence,
    and returns the p-values for each peak. If the largest significant
    peak is in the circadian range as specified by the args, it is
    rhythmic."""
    timer = plo.laptimer()
    print "Lomb-Scargle Periodogram... time: ",
    rhythmic_or_not = np.zeros(len(ls_data.T))
    pgram_data = np.zeros((300,len(ls_data.T)))
    circadian_peaks = np.zeros(len(ls_data.T))
    circadian_peak_periods = np.zeros(len(ls_data.T))
    # pgram
    for data_idx, d1 in enumerate(ls_data.T):
        # remove nans
        t1 = np.copy(times[~np.isnan(d1)])
        d1 = np.copy(d1[~np.isnan(d1)])
        pers, pgram, sig = blu.periodogram(t1, d1, period_low=1,
                        period_high=60, res=300)
        peak = np.argmax(pgram)

        if (pers[peak] >= circ_low and pers[peak] <= circ_high and
                                                    sig[peak] <=alpha):
            rhythmic_or_not[data_idx] = 1
            #
            circadian_peaks[data_idx]= pgram[peak]
            circadian_peak_periods[data_idx]= pers[peak]
        else:
            minpeak = np.argmax(pers>=circ_low)
            maxpeak = np.argmin(pers<circ_high)
            circadian_peaks[data_idx] = np.max(pgram[minpeak:maxpeak])
            circadian_peak_periods[data_idx] =\
                pers[minpeak:maxpeak][np.argmax(pgram[minpeak:maxpeak])]

        # return either normed or un-normed data
        pgram_data[:,data_idx] = pgram
    print str(np.round(timer(),1))+"s"
    return pers, pgram_data, circadian_peaks, circadian_peak_periods, rhythmic_or_not

def eigensmooth(times, detrended_data, ev_threshold=0.05, dim=40, min_ev=2):
    """
    Uses an eigendecomposition to keep only elements with >threshold of the
    data. Then it returns the denoised data.

    Notes: This should take the covariance matrix of the data using fwd-backward method. Then it
    eigendecomposes it, then it finds the biggest (2+) eigenvalues and returns only the
    components associated with them.
    For an intuitive explanation of how this works:
    http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    """
    timer = plo.laptimer()
    print "Eigendecomposition... time: ",

    #set up results - by default set to NaNs
    denoised_data = np.nan*np.ones(detrended_data.shape)
    eigenvalues_list = []

    # denoise each piece of data
    for data_idx, d0 in enumerate(detrended_data.T):
        # remove nans from times, d1. keep nans in d0
        t1 = times[~np.isnan(d0)]
        d1 = np.copy(d0[~np.isnan(d0)])

        # using spectrum to get the covariance matrix
        X = spectrum.linalg.corrmtx(d1, dim-1, method='autocorrelation')
        # the embedding matrix
        X = (1/np.sqrt(len(X.T)))*np.array(X)
        XT = np.transpose(X)
        C = XT.dot(X)

        # now eigendecompose
        evals, Q = np.linalg.eig(C)

        # find evals that matter, use a minimum of 2
        eval_goodness = np.max([2,
                        np.sum(evals/np.sum(evals) >= ev_threshold)])
        QT = np.transpose(Q)

        # and return the reconstruction
        P = QT.dot(XT)
        denoised = np.sum(P[:eval_goodness],0)

        # find alignment - for some reason the signal can be flipped.
        # this catches it
        # truncate the first 24h during alignment
        align, atype = alignment(d1, denoised, d=dim, dstart=96)

        # fix alignment if leading nans
        nanshift=0
        if np.isnan(d0[0]):
            nanshift=np.argmin(np.isnan(d0))

        # get the correctly-shaped denoised data
        denoised = denoised[align:align+len(d1)]*atype

        #align denoised data in matrix
        denoised_data[nanshift:len(denoised)+nanshift, data_idx] = denoised
        eigenvalues_list.append(evals/np.sum(evals))

    print str(np.round(timer(),1))+"s"
    return times, denoised_data, eigenvalues_list

def sinusoidal_fitting(times, data, rhythmic_or_not, fit_times=None,
                       forced_periods=None):
    """
    Takes detrended and denoised data and times, and fits a damped sinusoid to
    the data. Returns: times, sine_data, phase_data, periods, amplitudes,
    decay parameters, and pseudo rsq values.

    Has the option to set "full_times", the complete set of times the sinusoid
    should be fit for.

    If forced_periods are supplied, the sine fit will be within 1h of the forced
    period.
    """
    timer = plo.laptimer()
    print "Sinusoidal fitting... time: ",

    if fit_times is None:
        fit_times = np.copy(times)

    # these will be the outputs
    sine_data = np.nan*np.ones((len(fit_times), len(data[0])))
    phase_data = np.nan*np.ones((len(fit_times), len(data[0])))
    phases = np.nan*np.ones(len(data[0]))
    meaningful_phases = np.nan*np.ones(len(data[0]))
    periods = np.nan*np.ones(len(data[0]))
    amplitudes = np.nan*np.ones(len(data[0]))
    decays = np.nan*np.ones(len(data[0]))
    pseudo_r2s = np.nan*np.ones(len(data[0]))

    for idx,idata in enumerate(data.T):
            # copy data for editing
            cell = np.copy(idata)
            model = dsin.DecayingSinusoid(times, cell)
            # fit the data with the polynomial+sinusoid
            # note we are not using model averaging, we are
            # just taking the single best model by AICc

            # if no estimate is given just do the normal fitting
            # if an estimate is given, bound the period within two
            model.run()
            params = model.best_model.result.params

            if forced_periods is not None:
                # only use force if necessary
                if np.abs(params['period'].value-forced_periods[idx]) >1:
                    # force is necessary
                    model._estimate_parameters()
                    model._fit_models(period_force = forced_periods[idx])
                    model._calculate_averaged_parameters()
                    params = model.best_model.result.params
            # put the sine fit data in the sine_data matrix,
            # and same for phase
            # export all phases for heatmap, and other data for rhythmic cells only
            phases[idx] = (fit_times[0]*2*np.pi/params['period']+params['phase'])%(2*np.pi)

            if rhythmic_or_not[idx]==1:
                phase_data[:,idx] = (fit_times*2*np.pi/params['period']+params['phase'])%(2*np.pi)
                sine_data[:,idx] = dsin.sinusoid_component(params, fit_times)
            # summary stats
                periods[idx] = model.best_model.result.params['period']
                amplitudes[idx] = model.best_model.result.params['amplitude']
                decays[idx] = model.best_model.result.params['decay']
                pseudo_r2s[idx] = model.best_model._calc_r2()
                meaningful_phases[idx] = (fit_times[0]*2*np.pi/params['period']+params['phase'])%(2*np.pi)

    print str(np.round(timer(),1))+"s"
    return fit_times, sine_data, phase_data, phases, periods, amplitudes, decays, pseudo_r2s, meaningful_phases

def plot_result(cellidx, raw_times, raw_data, trendlines, detrended_times, detrended_data, eigenvalues, final_times, final_data, rhythmic_or_not, pers, pgram_data, sine_times, sine_data, r2s, output_folder, data_type):
    """
    Plotting utility to plot ALL data at once. A bit of a mess.
    """
    plo.PlotOptions(ticks='in')
    fig = plt.figure(figsize=(4,7))
    gs = gridspec.GridSpec(4,2)

    ax = plt.subplot(gs[0,0])
    bx = plt.subplot(gs[0,1])
    cx = plt.subplot(gs[1,0])
    dx = plt.subplot(gs[1,1])
    ex = plt.subplot(gs[2,0])
    fx = plt.subplot(gs[2,1])
    gx = plt.subplot(gs[3,0])
    hx = plt.subplot(gs[3,1])

    ax.plot(raw_times, raw_data[:,cellidx], label='raw', color='k')
    ax.set_ylabel('Raw Biolum (AU)')
    ax.set_xlim([0,np.max(raw_times)])

    bx.plot(raw_times, raw_data[:,cellidx], label='raw', color='k')
    valid = trendlines[:,cellidx]>0
    bx.plot(detrended_times[valid], trendlines[valid,cellidx],
                label='trend', color='gold')
    bx.set_xlim([0,np.max(raw_times)])

    cx.plot(detrended_times, detrended_data[:,cellidx],color='gold')
    cx.axhline(0,color='k', ls=':')
    cx.set_xlim([0,np.max(raw_times)])
    cx.set_ylabel('Detrended (AU)')

    dx.plot(eigenvalues[cellidx], marker='o',color='gold')
    dx.axhline(0.05, color='h', ls='--', label='cutoff')
    dx.set_ylabel('EVal (% of total)')

    ex.plot(final_times, final_data[:,cellidx],'h')
    ex.axhline(0,color='k', ls=':')
    ex.set_ylabel('Denoised (AU)')
    ex.set_xlim([0,np.max(raw_times)])

    if rhythmic_or_not[cellidx]==1:
        rhyth='RHY'
    else:
        rhyth='ARR'
    fx.plot(pers, pgram_data[:,cellidx])
    fx.set_ylim([0,1])
    fx.set_xlabel('Period')
    fx.set_ylabel('LS Pgram: '+rhyth)

    gx.plot(sine_times, sine_data[:,cellidx], 'f')
    gx.set_ylabel('Cosine, R$^2$= '+str(np.round(r2s[cellidx],2)))
    gx.set_xlim([0,np.max(raw_times)])

    hx.plot(detrended_times, detrended_data[:,cellidx],color='gold')
    hx.plot(final_times, final_data[:,cellidx],'h')
    hx.plot(sine_times, sine_data[:,cellidx], 'f')
    hx.set_xlim([0,np.max(raw_times)])
    hx.set_ylabel('Biolum (AU)')

    plt.tight_layout(**plo.layout_pad)
    plt.savefig(output_folder+data_type+'_cell'+str(cellidx)+'.png', dpi=300)
    plt.close()
    plt.clf()
