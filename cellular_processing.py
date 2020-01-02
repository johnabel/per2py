from __future__ import division

# imports
import numpy  as np
import scipy as sp
import pandas as pd
from matplotlib import gridspec
import matplotlib.pyplot as plt

# local functions to import
from LocalImports import PlotOptions as plo
from LocalImports import Bioluminescence as blu
from LocalImports import DecayingSinusoid as dsin
from LocalImports import CellularFunctions as cpf

#inputs nms pre
pull_from_imagej = True
input_folder = 'Demo/cellular/scn15_NMSWT_012218/' # edit this
input_extension = '.csv'

# do we want plots?
plotcount = 2
# for each dataset, this is formatted [descriptor, root filename, color ('Red', or 'Green')]
input_ij_files   = ['012218NMS_Green_Pre_Spots in tracks statistics',
                    '012218NMS_Red_Pre_Spots in tracks statistics',
                    '012218NMS_Green_TTX_Spots in tracks statistics',
                    '012218NMS_Red_TTX_Spots in tracks statistics',
                    '012218NMS_Green_Wash_Spots in tracks statistics',
                    '012218NMS_Red_Wash_Spots in tracks statistics'] # edit this

# list all the datasets
all_inputs=[]
for input_ij in input_ij_files:
    all_inputs.append(cpf.generate_filenames_dict(input_folder, input_ij, 
                                    pull_from_imagej, input_ij_extension=input_extension))

# process the data for every set of inputs
for files_dict in all_inputs:
    
    # assign all filenames to correct local variables
    data_type = files_dict['data_type']
    input_data = files_dict['input_data']
    input_folder = files_dict['input_folder']
    input_ij_extension = files_dict['input_ij_extension']
    input_ij_file = files_dict['input_ij_file']
    output_cosine = files_dict['output_cosine']
    output_cosine_params = files_dict['output_cosine_params']
    output_detrend = files_dict['output_detrend']
    output_zscore = files_dict['output_zscore']
    output_detrend_smooth = files_dict['output_detrend_smooth']
    output_detrend_smooth_xy = files_dict['output_detrend_smooth_xy']
    output_pgram = files_dict['output_pgram']
    output_phases = files_dict['output_phases']
    pull_from_imagej = files_dict['pull_from_imagej']
    raw_signal = files_dict['raw_signal']
    raw_xy = files_dict['raw_xy']    

    # does the actual processing of the data
    # I. IMPORT DATA
    # only perform this step if pull_from_imagej is set to True
    if pull_from_imagej:
        cpf.load_imagej_file(input_data, raw_signal, raw_xy)

    raw_times, raw_data, locations, header = cpf.import_data(raw_signal, raw_xy)

    # II. INTERPOLATE MISSING PARTS
    # truncate 0 h and interpolate
    interp_times, interp_data, locations = cpf.truncate_and_interpolate(
        raw_times, raw_data, locations, truncate_t=0)

    # III. DETREND USING HP Filter
    #(Export data for presentation of raw tracks with heatmap in Prism.)
    detrended_times, detrended_data, trendlines = cpf.hp_detrend(
                                        interp_times, interp_data)

    # IV. SMOOTHING USING EIGENDECOMPOSITION
    # eigendecomposition
    denoised_times, denoised_data, eigenvalues = cpf.eigensmooth(detrended_times,
        detrended_data, ev_threshold=0.05, dim=40)
    # TRUNCATE 12 INITIAL HOURS
    final_times, final_data, locations = cpf.truncate_and_interpolate(denoised_times,
                                    denoised_data, locations, truncate_t=12)

    # V. LS PERIODOGRAM TEST FOR RHYTHMICITY
    lspers, pgram_data, circadian_peaks, lspeak_periods, rhythmic_or_not = cpf.LS_pgram(final_times, final_data)

    # VI. GET A SINUSOIDAL FIT TO EACH CELL
    # use final_times, final_data
    # use forcing to ensure period within 1h of LS peak period
    sine_times, sine_data, phase_data, refphases, periods, amplitudes, decays, r2s, meaningful_phases =\
         cpf.sinusoidal_fitting(final_times, final_data, rhythmic_or_not, 
                               fit_times=raw_times, forced_periods=lspeak_periods)
    # get metrics
    circadian_metrics = np.vstack([rhythmic_or_not, circadian_peaks, refphases, periods, amplitudes,
                                   decays, r2s])

    # VII. SAVING ALL COMPONENTS
    timer = plo.laptimer()
    print "Saving data... time: ",

    # detrended
    cell_ids = header[~np.isnan(header)]
    output_array_det = np.nan*np.ones((len(detrended_times)+1, len(cell_ids)+2))
    output_array_det[1:,0] = detrended_times
    output_array_det[1:,1] = np.arange(len(detrended_times))
    output_array_det[0,2:] = refphases
    output_array_det[1:,2:] = detrended_data
    output_df = pd.DataFrame(data=output_array_det,
            columns = ['TimesH', 'Frame']+list(cell_ids))
    output_df.loc[0,'Frame']='RefPhase'
    output_df.to_csv(output_detrend, index=False)
    del output_df # clear it

    # detrended-denoised
    output_array = np.nan*np.ones((len(final_times)+1, len(cell_ids)+2))
    output_array[1:,0] = final_times
    output_array[1:,1] = np.arange(len(final_times))
    output_array[0,2:] = refphases
    output_array[1:,2:] = final_data
    output_df = pd.DataFrame(data=output_array,
            columns = ['TimesH', 'Frame']+list(cell_ids))
    output_df.loc[0,'Frame']='RefPhase'
    output_df.to_csv(output_detrend_smooth, index=False)
    del output_df # clear it
    
    # Z-Score
    output_array = np.nan*np.ones((len(final_times)+1, len(cell_ids)+2))
    output_array[1:,0] = final_times
    output_array[1:,1] = np.arange(len(final_times))
    output_array[1:,2:] = sp.stats.zscore(final_data, axis=0, ddof=0)
    output_df = pd.DataFrame(data=output_array,
            columns = ['TimesH', 'Frame']+list(cell_ids))
    output_df.loc[0,'Frame']='RefPhase'
    output_df.loc[0,list(cell_ids)]=refphases
    output_df.to_csv(output_zscore, index=False)
    del output_df # clear it

    # LS Pgram
    output_array = np.nan*np.ones((len(lspers), len(pgram_data[0,:])+1))
    output_array[:,0] = lspers
    output_array[:,1:] = pgram_data
    output_df = pd.DataFrame(data=output_array,
            columns = ['LSPeriod']+list(cell_ids))
    output_df.to_csv(output_pgram, index=False)
    del output_df # clear it

    #sinusoids
    output_array = np.nan*np.ones((len(sine_times), len(cell_ids)+2))
    output_array[:,0] = sine_times
    output_array[:,1] = np.arange(len(sine_times))
    output_array[:,2:] = sine_data
    output_df = pd.DataFrame(data=output_array,
            columns = ['TimesH', 'Frame']+list(cell_ids))
    output_df.to_csv(output_cosine, index=False)
    del output_df

    #phases
    output_array = np.nan*np.ones((len(sine_times), len(cell_ids)+2))
    output_array[:,0] = sine_times
    output_array[:,1] = np.arange(len(sine_times))
    output_array[:,2:] = phase_data
    output_df = pd.DataFrame(data=output_array,
            columns = ['TimesH', 'Frame']+list(cell_ids))
    output_df.to_csv(output_phases, index=False)
    del output_df

    # sinusoid parameters and XY locations
    # this gets the locations for each cell by just giving their mean
    # location and ignoring the empty values. this is a fine approximation.
    locs_fixed = np.zeros([2,len(cell_ids)])
    for idx in range(len(cell_ids)):
        locs_fixed[0, idx] = np.nanmean(locations[:,idx*2])
        locs_fixed[1, idx] = np.nanmean(locations[:,idx*2+1]) 
    output_array = np.nan*np.ones((9, len(cell_ids)))
    output_array= np.concatenate((circadian_metrics,locs_fixed), axis=0)
    output_array[2,:] *= 360/2/np.pi #transform phase into 360-degree circular format
    output_df = pd.DataFrame(data=output_array,
            columns = list(cell_ids), index=['Rhythmic','CircPeak','Phase','Period','Amplitude',
                                            'Decay','Rsq', 'X', 'Y'])
    output_df.T.to_csv(output_cosine_params, index=True)
    del output_df # clear it
    print str(np.round(timer(),1))+"s"
    
    print "Generating and saving plots: ",
    cellidxs=np.random.randint(len(cell_ids),size=plotcount)
    for cellidx in cellidxs:
        # truly awful syntax
        cpf.plot_result(cellidx, raw_times, raw_data, trendlines, 
                    detrended_times, detrended_data, eigenvalues, 
                    final_times, final_data, rhythmic_or_not, 
                    lspers, pgram_data, sine_times, sine_data, r2s,
                    output_folder, data_type)
    print str(np.round(timer(),1))+"s"

    print "All data saved. Run terminated successfully for "+data_type+'.\n'