"""
Created on Mon April 1 2019

@author: john abel

Module set up to perform analysis for whole-body circadian recordings
by Yongli Shan.
"""
from __future__ import division

import numpy  as np
import pandas as pd
from scipy import signal, interpolate, optimize, sparse, stats
from scipy.sparse import dia_matrix, eye as speye
from scipy.sparse.linalg import spsolve
import spectrum
import matplotlib.pyplot as plt

import collections

from . import PlotOptions as plo
from . import DecayingSinusoid as ds
from . import Bioluminescence as blu


class WholeBodyRecording(object):
    """
    Class to analyze time series data for whole body Per2iLuc recordings.
    """

    def __init__(self, directory, red_file, green_file, imaging_start, imaging_interval, name=None):
        """
        Parameters
        ----------
        directory : str        
            The path to the folder that contains the data files.
        red_data : np.ndarray
            Should contain a single PER2iLuc time series for analysis.
        green_data : np.ndarray
            Should contain a single PER2iLuc time series for analysis.
        imaging_start : str
            Date-time start of the imaging, in format '2019-03-12 17:08:01'.
        imaging_interval : str
            Timespan of each image, in format 'XYZ min'.
        name : str, optional
            A name for the dataset if desired.
        """

        if name is not None:
            self.name = name
        self.directory = directory+'/'


        red_f = np.genfromtxt(self.directory+red_file)
        green_f = np.genfromtxt(self.directory+green_file)
        imaging_times = red_f[1:,0]
        assert len(imaging_times)==len(green_f[1:,0]), "Imaging files of unequal length."

        # do we need to remove outliers? or can we just use the LSPgram to get it.....
        ry = red_f[1:,1]
        gy = green_f[1:,1]
        xi = pd.date_range(imaging_start, periods=len(ry), freq=imaging_interval)

        imaging = {}
        imaging['red'] = ry
        imaging['green'] = gy
        imaging['t_imaging'] = xi
        self.imaging = imaging
        self.imaging_start = imaging_start
        self.imaging_interval = imaging_interval
        self.corrmat = {}
    
    def excise_imaging_times(self, intervals, t_pre='t_imaging', red_pre='red', green_pre='green', t_post='t_excised', red_post='red_excised', green_post='green_excised', cooldown_ext=5):
        """Cuts times where the signal is interfered with due to light in the animal's cage. By default, it operates on and modifies t_imaging, red, green, and creates t_excised, red_excised, and green_excised in the self.imaging dict.
        
        Parameters
        ----------
        intervals : list of lists of length 2
            A list of lists, where each list contains the times and dates between which the data should be excised. For example, [['2019-03-01 07:52:00', '2019-03-01 08:24:00']].
        t_pre : str, optional
            Key for timeseries in the self.imaging dictionary upon which the excision operates, by default 't_imaging'
        red_pre : str, optional
            Key for red series in the self.imaging dictionary upon which the excision operates, by default 'red'
        green_pre : str, optional
            Key for green series in the self.imaging dictionary upon which the excision operates, by default 'green'
        t_post : str, optional
            Key for placing the timeseries corresponding to the excised data into the self.imaging dict, by default 't_excised'
        red_post : str, optional
            Key for placing the excised red data into the self.imaging dict, by default 'red_excised'
        green_post : str, optional
            Key for placing the excised green data into the self.imaging dict, by default 'green_excised'
        cooldown_ext : int, optional
            Time, in minutes, where the data should be excised following the end of a light epoch due to artifact, by default 5
        """        

        # if pre is not set, use the already-truncated data
        # if already-truncated data does not exist, take the raw
        if t_pre not in self.imaging.keys():
            self.imaging[t_pre] = self.imaging['t_imaging']

        if red_pre not in self.imaging.keys():
            self.imaging[red_pre] = self.imaging['red']

        if green_pre not in self.imaging.keys():
            self.imaging[green_pre] = self.imaging['green']
        
        # get data for editing
        xr = self.imaging[t_pre]
        y1r = self.imaging[red_pre]
        y2r = self.imaging[green_pre]
        
        # cut each pulse
        for pulse in intervals:
            lightson  = pd.to_datetime(pulse[0])
            lightsoff = pd.to_datetime(pulse[1])
            try:
                idx = np.where(xr>=lightson)[0][0]
                idx2 = np.where(xr>=lightsoff)[0][0]+cooldown_ext
                y1r = np.hstack([y1r[:idx], y1r[idx2:]])
                y2r = np.hstack([y2r[:idx], y2r[idx2:]])
                xr = xr[:idx].append(xr[idx2:])
            except IndexError:
                # if there is no data that late
                pass

        # by default drop it in these
        if t_post is None:
            t_post = 't_excised'
        if red_post is None:
            red_post = 'red_excised'
        if green_post is None:
            green_post is 'green_excised'

        self.imaging[red_post] = y1r
        self.imaging[green_post] = y2r
        self.imaging[t_post] = xr

    def import_temperature_humidity(self, filename, start_time='imaging', droplines=18):
        """Imports the temperature and humidity file and links to existing data. Creates self.TH, a dictionary of the temperature and humidity recordings.
        
        Parameters
        ----------
        filename : str
            Name of the temperature and humidity file, within self.directory.
        start_time : str, optional
            The time at which the temperature and humidity recordings start, in format '2019-03-12 17:08:01', by default, 'imaging', the time at which the imaging also starts.
        droplines : int, optional
            Number of lines to discard at the start of the file, by default 18 (which works for the data as in the Demo directory).
        """

        # load the files
        TH_pd = pd.read_excel(self.directory+filename, usecols=[0,1,2,3],
                              names=['index', 'date', 'temp', 'humidity'])
        TH_pd = TH_pd.drop(range(droplines))
        # if set up as example file, then this is 
        # the correct number of rows to drop

        # collect the x values for temp and humidity
        xth = pd.DatetimeIndex(pd.to_datetime(TH_pd['date'], yearfirst=True))

        # only start where imaging starts
        if start_time=='imaging':
            imgstart  = pd.to_datetime(self.imaging_start)
        elif start_time==None:
            imgstart = pd.to_datetime('2010-01-01 00:00:00')
        else:
            try:
                imgstart = pd.to_datetime(start_time)
            except:
                imgstart = pd.to_datetime('2010-01-01 00:00:00')
                print("Date format not understood, leaving all times.")
        idx = np.where(xth>=imgstart)[0][0]
        xthr = xth[idx:]
        tempr = np.array(TH_pd['temp'], dtype=np.float)[idx:]
        humr = np.array(TH_pd['humidity'], dtype=np.float)[idx:]

        interval = xthr[1]-xthr[0]
        offset = xthr[0]-imgstart

        TH = {}
        TH['temp'] = tempr
        TH['hum'] = humr
        TH['t'] = xthr
        TH['interval_h'] = interval.total_seconds()/60/60
        TH['offset_h'] = offset.total_seconds()/60/60
        self.TH = TH

    def import_actogram(self, filename, start_time='imaging', actogram_interval=1):
        """Imports the actogram file and attaches it to the object at self.actogram.

        Parameters
        ----------
        filename : str
            Name of the temperature and humidity file, within self.directory.
        start_time : str, optional
            The time at which the actogram start, in format '2019-03-12 17:08:01'. Defulats to 'imaging', the time at which the imaging also starts.
        actogram_interval : float, optional, defaults to 1
            The interval of the actogram recording in minutes.
        """
        act_pd = pd.read_csv(self.directory+filename, header=None)
        total_days = act_pd.shape[1]-1
        act_start = pd.to_datetime(act_pd[1][1][3:]+' 00:00:00', dayfirst=True)

        # only start where imaging starts
        if start_time=='imaging':
            imgstart  = pd.to_datetime(self.imaging_start)
        elif start_time==None:
            imgstart = pd.to_datetime('2010-01-01 00:00:00')
        else:
            try:
                imgstart = pd.to_datetime(start_time)
            except:
                imgstart = pd.to_datetime('2010-01-01 00:00:00')
                print("Date format not understood, leaving all times.")

        # assemble all the columns
        intervals = int(60/actogram_interval*24)*total_days
        xa = pd.date_range(act_start, periods=intervals, freq=str(actogram_interval)+' min')
        activ = np.array(
                    act_pd.iloc[np.arange(8,8+(60/actogram_interval*24)),
                    np.arange(1,1+total_days)], dtype=np.float
                             ).flatten('F')
        idx = np.where(xa>=imgstart)[0][0]
        xar = xa[idx:]
        actr = activ[idx:]
        offset = xar[0] - imgstart

        # remove nans
        xar = xar[~np.isnan(actr)]
        actr = actr[~np.isnan(actr)]

        interval_h = actogram_interval/60
        offset_h = offset.total_seconds()/60/60

        activity = {}
        activity['t'] = xar
        activity['activity'] = actr
        activity['interval_h'] = interval_h
        activity['offset_h'] = offset_h
        self.activity = activity


    def process_imaging_data(self, tname, redname, greenname, lsperiod_fit=False):
        """Performs the analysis of the PER2iLuc data. Consists of a Hodrick-Prescott detrend, a Lomb-Scargle periodogram test for rhythmicity, a smoothing via eigendecomposition, fitting to a cosine function, and returning relevant results in the self.imaging, self.periodogram, and self.sinusoid dictionaries. Importantly, this results in a timeseries with the ending '_UT' which denotes that it is in universal time (i.e., time since the start of the imaging).
        
        Parameters
        ----------
        tname : str
            Key for the timeseries in the self.imaging dict to be used.
        redname : str
            Key for the red series in the self.imaging dict to be used.
        greenname : str
            Key for the green series in the self.imaging dict to be used.
        lsperiod_fit : bool, optional
            If True, bounds the resulting sinusoid to have a period within 1h of LSPgram. By default False
        """
        x = self.imaging[tname]
        red = self.imaging[redname]
        green = self.imaging[greenname]

        # convert interval into float
        timediffs = [x[i+1]-x[i] for i in range(len(x)-1)]
        timediffs_h = [td.total_seconds()/60/60 for td in timediffs]
        times = np.hstack([0,np.cumsum(timediffs_h)])
        # save times that relate to other measurements
        self.imaging[tname+'_UT']  = times

        hpt, hp_red, hp_redb = hp_detrend(times, red)
        hpt, hp_green, hp_greenb = hp_detrend(times, green)

        self.imaging[redname+'_hp'] = hp_red
        self.imaging[greenname+'_hp'] = hp_green
        self.imaging[redname+'_hpb'] = hp_redb
        self.imaging[greenname+'_hpb'] = hp_greenb
        
        # Use LSPgram to confirm rhythmic or not
        red_pgram = circadian_LSPgram(hpt, hp_red, circ_low=18, 
                                      circ_high=30, alpha=0.05)
        green_pgram = circadian_LSPgram(hpt, hp_green, circ_low=18, 
                                      circ_high=30, alpha=0.05)
        try:
            self.periodogram['red'] = red_pgram
            self.periodogram['green'] = green_pgram
        except AttributeError:
            self.periodogram = {}
            self.periodogram['red'] = red_pgram
            self.periodogram['green'] = green_pgram            

        # now the eigensmooth
        et, ered, _ = eigensmooth(hpt, hp_red)
        et, egreen, _ = eigensmooth(hpt, hp_green)

        self.imaging[redname+'_es'] = ered
        self.imaging[greenname+'_es'] = egreen

        # fit the data with the polynomial+sinusoid
        # note we are not using model averaging, we are
        # just taking the single best model by AICc
        rmod = ds.DecayingSinusoid(et, ered)
        gmod = ds.DecayingSinusoid(et, egreen)

        # if no estimate is given just do the normal fitting
        # if an estimate is given, bound the period within 1h
        rmod.run()
        rparams = rmod.best_model.result.params
        gmod.run()
        gparams = gmod.best_model.result.params

        if lsperiod_fit:
            # only use force if necessary
            if np.abs(rparams['period'].value-red_pgram['period']) >1:
                # force is necessary
                rmod._estimate_parameters()
                rmod._fit_models(period_force=red_pgram['period'])
                rmod._calculate_averaged_parameters()
                rparams = rmod.best_model.result.params
            if np.abs(gparams['period'].value-green_pgram['period']) >1:
                # force is necessary
                gmod._estimate_parameters()
                gmod._fit_models(period_force=green_pgram['period'])
                gmod._calculate_averaged_parameters()
                gparams = gmod.best_model.result.params
        # put the sine fit data in the sine_data matrix,
        # and same for phase
        # export all phases for heatmap, and other data for rhythmic cells only
        # rphases = (et*2*np.pi/rparams['period']+rparams['phase'])%(2*np.pi)
        # gphases = (et*2*np.pi/gparams['period']+gparams['phase'])%(2*np.pi)

        try:
            self.sinusoids
        except AttributeError:
            self.sinusoids = {}

        if red_pgram['rhythmic']:
            red_sin = {}
            red_sin['ts'] = et
            red_sin['phase_data'] = (et*2*np.pi/rparams['period']+rparams['phase'])%(2*np.pi)
            red_sin['sine_data'] = ds.sinusoid_component(rparams, et)
            # summary stats
            red_sin['pseudo_r2'] = rmod.best_model._calc_r2()
            red_sin['params'] = rparams
            self.sinusoids['red'] = red_sin
        if green_pgram['rhythmic']:
            green_sin = {}
            green_sin['ts'] = et
            green_sin['phase_data'] = (et*2*np.pi/gparams['period']+gparams['phase'])%(2*np.pi)
            green_sin['sine_data'] = ds.sinusoid_component(gparams, et)
            # summary stats
            green_sin['pseudo_r2'] = gmod.best_model._calc_r2()
            green_sin['params'] = gparams
            self.sinusoids['green'] = green_sin


    def process_temp_hum_data(self, tname='t', tempname='temp', humname='hum', lsperiod_fit=False):
        """Performs the analysis of the temperature and humidity data. May need adjusting depending on file format. Consists of a Hodrick-Prescott detrend, a Lomb-Scargle periodogram test for rhythmicity, a smoothing via eigendecomposition, fitting to a cosine function, and returning relevant results in the self.TH, self.periodogram, and self.sinusoid dictionaries. Importantly, this results in a timeseries with the ending '_UT' which denotes that it is in universal time (i.e., time since the start of the imaging).
        
        Parameters
        ----------
        tname : str, optional
            Key for the timeseries in the self.imaging dict to be used, by default 't'
        tempname : str
            Key for the temperature series in the self.imaging dict to be used, by default 'temp'
        humname : str
            Key for the humidity series in the self.imaging dict to be used, by default 'hum'
        lsperiod_fit : bool, optional
            If True, bounds the resulting sinusoid to have a period within 1h of LSPgram. By default False
        """
        x = self.TH[tname]
        temp = self.TH[tempname]
        hum = self.TH[humname]

        # convert interval into float
        interval_float = self.TH['interval_h']
        times = np.arange(len(x))*interval_float+self.TH['offset_h'] # I think
        self.TH[tname+'_UT'] = times

        hpt, hp_temp, hp_tempb = hp_detrend(times, temp)
        hpt, hp_hum, hp_humb = hp_detrend(times, hum)

        self.TH[tempname+'_hp'] = hp_temp
        self.TH[humname+'_hp'] = hp_hum
        self.TH[tempname+'_hpb'] = hp_tempb
        self.TH[humname+'_hpb'] = hp_humb

        
        # Use LSPgram to confirm rhythmic or not
        temp_pgram = circadian_LSPgram(hpt, hp_temp, circ_low=18, 
                                      circ_high=30, alpha=0.05)
        hum_pgram = circadian_LSPgram(hpt, hp_hum, circ_low=18, 
                                      circ_high=30, alpha=0.05)
        try:
            self.periodogram['temp'] = temp_pgram
            self.periodogram['hum'] = hum_pgram
        except AttributeError:
            self.periodogram = {}
            self.periodogram['temp'] = temp_pgram
            self.periodogram['hum'] = hum_pgram            

        # now the eigensmooth
        et, etemp, _ = eigensmooth(hpt, hp_temp)
        et, ehum, _ = eigensmooth(hpt, hp_hum)

        self.TH[tempname+'_es'] = etemp
        self.TH[humname+'_es'] = ehum

        # fit the data with the polynomial+sinusoid
        # note we are not using model averaging, we are
        # just taking the single best model by AICc
        rmod = ds.DecayingSinusoid(et, etemp)
        gmod = ds.DecayingSinusoid(et, ehum)

        # if no estimate is given just do the normal fitting
        # if an estimate is given, bound the period within 1h
        rmod.run()
        rparams = rmod.best_model.result.params
        gmod.run()
        gparams = rmod.best_model.result.params

        if lsperiod_fit:
            # only use force if necessary
            if np.abs(rparams['period'].value-temp_pgram['period']) >1:
                # force is necessary
                rmod._estimate_parameters()
                rmod._fit_models(period_force=temp_pgram['period'])
                rmod._calculate_averaged_parameters()
                rparams = rmod.best_model.result.params
            if np.abs(gparams['period'].value-hum_pgram['period']) >1:
                # force is necessary
                gmod._estimate_parameters()
                gmod._fit_models(period_force=hum_pgram['period'])
                gmod._calculate_averaged_parameters()
                gparams = gmod.best_model.result.params
        # put the sine fit data in the sine_data matrix,
        # and same for phase
        # export all phases for heatmap, and other data for rhythmic cells only
        # rphases = (et*2*np.pi/rparams['period']+rparams['phase'])%(2*np.pi)
        # gphases = (et*2*np.pi/gparams['period']+gparams['phase'])%(2*np.pi)

        try:
            self.sinusoids
        except AttributeError:
            self.sinusoids = {}

        if temp_pgram['rhythmic']:
            temp_sin = {}
            temp_sin['ts'] = et
            temp_sin['phase_data'] = (et*2*np.pi/rparams['period']+rparams['phase'])%(2*np.pi)
            temp_sin['sine_data'] = ds.sinusoid_component(rparams, et)
            # summary stats
            temp_sin['pseudo_r2'] = rmod.best_model._calc_r2()
            temp_sin['params'] = rparams
            self.sinusoids['temp'] = temp_sin
        if hum_pgram['rhythmic']:
            hum_sin = {}
            hum_sin['ts'] = et
            hum_sin['phase_data'] = (et*2*np.pi/gparams['period']+gparams['phase'])%(2*np.pi)
            hum_sin['sine_data'] = ds.sinusoid_component(gparams, et)
            # summary stats
            hum_sin['pseudo_r2'] = gmod.best_model._calc_r2()
            hum_sin['params'] = gparams
            self.sinusoids['hum'] = hum_sin

    def process_activity_data(self, tname='t', actname='activity', lsperiod_fit=False, binsize=None):
        """Performs the analysis of the activity data.
        
        Parameters
        ----------
        tname : str, optional
            Key for the time series in the self.activity dict to be used,by default 't'
        actname : str, optional
            Key for the activity series in the self.activity dict to be used, by default 'activity'
        lsperiod_fit : bool, optional
            If True, bounds the resulting sinusoid to have a period within 1h of LSPgram. By default False
        binsize : int, optional
            Number of bins in each epoch, by default None. The circadian phase error for 15-min bins is at most 0.065 rad.
        """        
        x = self.activity[tname]
        act = self.activity[actname]

        # convert interval into float
        timediffs = [x[i+1]-x[i] for i in range(len(x)-1)]
        timediffs_h = [td.total_seconds()/60/60 for td in timediffs]
        times = np.hstack([0,np.cumsum(timediffs_h)]) +self.activity['offset_h']

        # bin the data if binning
        if type(binsize) is int:
            new_times = times[::binsize]
            bins = new_times
            digitized = np.digitize(times, bins)
            new_act = np.array([act[digitized == i].sum()
                            for i in range(1, len(bins))])
            # re-assign
            act = new_act
            times = new_times[1:] # time at *end* of bin
            self.activity[actname+'_b'] = act



        self.activity[tname+'_UT'] = times

        # get and remove baseline
        actb = np.mean(act)
        act_zero = act-actb

        self.activity[actname+'_zero'] = act_zero
        self.activity[actname+'_mean'] = actb
        
        # Use LSPgram to confirm rhythmic or not
        act_pgram = circadian_LSPgram(times, act_zero, circ_low=18, 
                                      circ_high=30, alpha=0.05)
        try:
            self.periodogram['act'] = act_pgram
        except AttributeError:
            self.periodogram = {}
            self.periodogram['act'] = act_pgram

        # now the eigensmooth
        et, eact, _ = eigensmooth(times, act_zero)
        self.activity[actname+'_es'] = eact

        # fit the data with the polynomial+sinusoid
        # note we are not using model averaging, we are
        # just taking the single best model by AICc
        rmod = ds.DecayingSinusoid(et, eact)

        # if no estimate is given just do the normal fitting
        # if an estimate is given, bound the period within 1h
        rmod.run()
        rparams = rmod.best_model.result.params

        if lsperiod_fit:
            # only use force if necessary
            if np.abs(rparams['period'].value-act_pgram['period']) >1:
                # force is necessary
                rmod._estimate_parameters()
                rmod._fit_models(period_force=act_pgram['period'])
                rmod._calculate_averaged_parameters()
                rparams = rmod.best_model.result.params

        # put the sine fit data in the sine_data matrix,
        # and same for phase
        # export all phases for heatmap, and other data for rhythmic cells only
        # rphases = (et*2*np.pi/rparams['period']+rparams['phase'])%(2*np.pi)

        try:
            self.sinusoids
        except AttributeError:
            self.sinusoids = {}

        if act_pgram['rhythmic']:
            act_sin = {}
            act_sin['ts'] = et
            act_sin['phase_data'] = (et*2*np.pi/rparams['period']+rparams['phase'])%(2*np.pi)
            act_sin['sine_data'] = ds.sinusoid_component(rparams, et)
            # summary stats
            act_sin['pseudo_r2'] = rmod.best_model._calc_r2()
            act_sin['params'] = rparams
            self.sinusoids['act'] = act_sin

    def continuous_wavelet_transform(self, dtype='es', shortestperiod=16, longestperiod=32, nvoice=512, be=5):
        """
        Gives CWT, phase, and period for all data marked with the data label.
        Tries to find r/g/t/h/a and applies CWT to each.
        """

        self.cwt = {}
        # imaging data
        for datatype in ['green_excised_'+dtype, 'red_excised_'+dtype]:
            x = self.imaging['t_excised_UT']
            try:
                y = self.imaging[datatype]
                cwt = blu.continuous_wavelet_transform(x,y,shortestperiod=shortestperiod, longestperiod=longestperiod, nvoice=nvoice, be=be)
                self.cwt[datatype] = cwt
            except:
                print("CWT not performed for "+datatype)
                pass

        # t/h data
        for datatype in ['temp_'+dtype, 'hum_'+dtype]:
            x = self.TH['x_UT']
            try:
                y = self.TH[datatype]
                cwt = blu.continuous_wavelet_transform(x,y,shortestperiod=shortestperiod, longestperiod=longestperiod, nvoice=nvoice, be=be)
                self.cwt[datatype] = cwt
            except:
                print("CWT not performed for "+datatype)
                pass

        # activity data
        try:
            # throws an optimization error
            x = self.activity['x_UT']
            y = self.activity['activity_es']
            cwt = blu.continuous_wavelet_transform(x,y,shortestperiod=shortestperiod, longestperiod=longestperiod, nvoice=nvoice, be=be)
            self.cwt['activity_es'] = cwt
        except:
            pass

    def plot_cwt_simple(self, dname='temp_es', name='', ax=None, colorbar=True, legend=True):
        """A simple plot of the continuous wavelet transform of a dataseries.
        
        Parameters
        ----------
        dname : str, optional
            [description], by default 'temp_es'
        name : str, optional
            [description], by default ''
        ax : [type], optional
            [description], by default None
        colorbar : bool, optional
            [description], by default True
        legend : bool, optional
            [description], by default True
        
        Returns
        -------
        [type]
            [description]
        """        

        if ax is None:
            ax = plt.subplot()

        if not hasattr(self, 'cwt'):
            self.continuous_wavelet_transform()

        cb = ax.pcolormesh(self.cwt[dname]['t'], self.cwt[dname]['tau'], 
                      self.cwt[dname]['cwt_scale'], cmap='jet')
        ax.set_xlabel('Time')
        ax.set_ylabel('Period')
        ax.plot(self.cwt[dname]['t'], self.cwt[dname]['period'], c='k',
                label='CWT Tau '+name)
        if colorbar:
            plt.colorbar(cb)
        if legend:
            plt.legend()
        return ax
    
    def correlate_signals(self, signals, dtypes, metric='pearsonr', max_dist=0.25, return_downsampled_trajectories=False, tmin=None, tmax=np.infty, name='correlation'):
        """
        Correlates signals using different metrics of correlation. Uses
        whichever signal is shorter as the reference, and minimizing the 
        time-distance between samples. Currently only uses Pearson's r.

        Arguments
        ----------
        signals : list
            The signals to correlate. Options: 'activity', 'biolum', 'TH'
        dtypes : list of str
            The type of data the analysis is performed on. Must be the same length as signals.
        metric : 'pearsonr'
            Statistic used for correlation. Currently only pearson's r is implemented. Others may be added simply here.
        max_dist : float, defaults to 0.25
            Maximum time distance to say timepoints are identical, in h. 
        """

        ts = [] # timeseries
        ds = [] # dataseries
        ss = [] # sorted signals

        for si, sig in enumerate(signals):
            if sig=='activity':
                ti = self.activity['x_UT']
                ts.append(ti[np.logical_and(ti>tmin, ti<tmax)])
                ds.append(self.activity['activity'+dtypes[si]][np.logical_and(ti>tmin, ti<tmax)])
                ss.append('activity'+dtypes[si])
            if sig=='biolum':
                ti = self.imaging['t_excised_UT']
                gi = self.imaging['green_excised'+dtypes[si]]
                ri = self.imaging['red_excised'+dtypes[si]]
                assert len(ti)==len(gi)==len(ri), "Time mismatch in biolum data."
                ts.append(ti[np.logical_and(ti>tmin, ti<tmax)])
                ts.append(ti[np.logical_and(ti>tmin, ti<tmax)])
                ds.append(gi[np.logical_and(ti>tmin, ti<tmax)])
                ds.append(ri[np.logical_and(ti>tmin, ti<tmax)])
                ss.append('biolum-green')
                ss.append('biolum-red')
            if sig=='TH':
                ti = self.TH['x_UT']
                ts.append(ti[np.logical_and(ti>tmin, ti<tmax)])
                ts.append(ti[np.logical_and(ti>tmin, ti<tmax)])
                ds.append(self.TH['temp'+dtypes[si]][np.logical_and(ti>tmin, ti<tmax)])
                ds.append(self.TH['hum'+dtypes[si]][np.logical_and(ti>tmin, ti<tmax)])
                ss.append('temp')
                ss.append('hum')

        corrmat = np.zeros([len(ds), len(ds)])
        pmat = np.ones([len(ds), len(ds)])

        for d1i in np.arange(len(ds)):
            for d2i in np.arange(len(ds)):
                if d1i <= d2i:

                    # get which data we are correlating
                    t1 = ts[d1i]
                    t2 = ts[d2i]
                    d1 = ds[d1i]
                    d2 = ds[d2i]
                    t1t2 = [t1, t2]
                    d1d2 = [d1, d2]

                    samelist = False
                    if len(t1)==len(t2):
                        if all(t1==t2):
                            samelist = True
                    
                    if samelist==True:
                        if metric=='pearsonr':
                            corrmat[d1i, d2i], pmat[d1i, d2i] = stats.pearsonr(d1, d2)
                        else:
                            print "Method of correlation not implemented."

                    else:
                        # reference is the shorter one, relative is the longer one
                        ref_ind = np.argmin([len(t1), len(t2)])
                        rel_ind = np.argmax([len(t1), len(t2)])

                        ref_t = t1t2[ref_ind]
                        ref_d = d1d2[ref_ind]
                        rel_t = t1t2[rel_ind]
                        rel_d = d1d2[rel_ind]
                        
                        rel_t_match = []
                        ref_t_match = []
                        rel_d_match = []
                        ref_d_match = []
                        for ref_idx,ti in enumerate(ref_t):
                            match_idx = np.argmin((rel_t-ti)**2)
                            if np.abs(rel_t[match_idx]-ti) < max_dist:
                                ref_t_match.append(ti)
                                ref_d_match.append(ref_d[ref_idx])
                                rel_t_match.append(rel_t[match_idx])
                                rel_d_match.append(rel_d[match_idx])

                        # do the correlations
                        if metric=='pearsonr':
                            corrmat[d1i, d2i], pmat[d1i, d2i] = stats.pearsonr(ref_d_match, rel_d_match)
                        else:
                            print "Method of correlation not implemented." 

        self.corrmat[name] = {'ps': pmat, 'corr': corrmat, 'names':ss}    




# functions
def hp_detrend(x, y, est_period=24., ret="both", a=0.05):
    """ Detrend the data using a hodrick-prescott filter. If ret ==
    "mean", return the detrended mean of the oscillation. Estimated
    period and 'a' parameter are responsible for determining the optimum
    smoothing parameter """

    x = np.asarray(x)
    y = np.asarray(y)

    # yt, index = timeseries_boundary(y, opt_b='mir', bdetrend=False)

    # As recommended by Ravn, Uhlig 2004, a calculated empirically
    num_periods = (x.max() - x.min())/est_period
    points_per_period = len(x)/num_periods
    w = a*points_per_period**4


    y_mean = hpfilter(y, w)
    y_detrended = y - y_mean

    if ret == "detrended": return x, y_detrended
    elif ret == "mean": return x, y_mean
    elif ret == "both": return x, y_detrended, y_mean


def hpfilter(X, lamb):
    """ Code to implement a Hodrick-Prescott with smoothing parameter
    lambda. Code taken from statsmodels python package (easier than
    importing/installing, https://github.com/statsmodels/statsmodels """

    X = np.asarray(X, float)
    if X.ndim > 1:
        X = X.squeeze()
    nobs = len(X)
    I = speye(nobs,nobs)
    offsets = np.array([0,1,2])
    data = np.repeat([[1.],[-2.],[1.]], nobs, axis=1)
    K = dia_matrix((data, offsets), shape=(nobs-2,nobs))

    trend = spsolve(I+lamb*K.T.dot(K), X, use_umfpack=True)
    return trend


def eigensmooth(times, data, ev_threshold=0.05, dim=600, min_ev=2):
    """
    Uses an eigendecomposition to keep only elements with >threshold of the
    data. Then it returns the denoised data.

    Notes: This should take the covariance matrix of the data using fwd-backward method. Then it
    eigendecomposes it, then it finds the biggest (2+) eigenvalues and returns only the
    components associated with them.
    For an intuitive explanation of how this works:
    http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    """
    # remove nans from times, d1. keep nans in d0
    t1 = times[~np.isnan(data)]
    d1 = np.copy(data[~np.isnan(data)])

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
    if np.isnan(data[0]):
        nanshift=np.argmin(np.isnan(data))
    else: nanshift = 0

    # get the correctly-shaped denoised data
    denoised = denoised[nanshift+align:nanshift+align+len(d1)]*atype

    return times, denoised, evals


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


def butterworth_lowpass(x, y, cutoff_period=8., order=10):
    """ Filter the data with a lowpass filter, removing noise with a
    critical frequency corresponding to the number of hours specified by
    cutoff_period. Assumes a period of 24h, with data in x in the units
    of hours. """

    x = np.asarray(x)
    y = np.asarray(y)
    nyquist = (x[1] - x[0])/2.
    cutoff_freq = 1/((cutoff_period/(x.max() - x.min()))*len(x))

    b, a = signal.butter(order, cutoff_freq/nyquist)
    y_filt = signal.filtfilt(b, a, y)

    return x, y_filt

def circadian_LSPgram(times, data, circ_low=18, circ_high=30, alpha=0.05):
    """Calculates a LS periodogram for each data sequence,
    and returns the p-values for each peak. If the largest significant
    peak is in the circadian range as specified by the args, it is
    rhythmic."""

    t1 = np.copy(times[~np.isnan(data)])
    d1 = np.copy(data[~np.isnan(data)])
    pers, pgram, sig = blu.periodogram(t1, d1, period_low=1,
                    period_high=60, res=300)
    peak = np.argmax(pgram)

    if (pers[peak] >= circ_low and pers[peak] <= circ_high and
                                                sig[peak] <=alpha):
        rhythmic_or_not = 1
        #
        circadian_peak = pgram[peak]
        circadian_peak_period = pers[peak]
    else:
        minpeak = np.argmax(pers>=circ_low)
        maxpeak = np.argmin(pers<circ_high)
        circadian_peak = np.max(pgram[minpeak:maxpeak])
        circadian_peak_period =\
            pers[minpeak:maxpeak][np.argmax(pgram[minpeak:maxpeak])]

    periodogram = {}
    periodogram['pers'] = pers
    periodogram['pgram'] = pgram
    periodogram['circadian_peak'] = circadian_peak
    periodogram['rhythmic'] = rhythmic_or_not
    periodogram['period'] = circadian_peak_period
    return periodogram

def periodogram(x, y, period_low=1, period_high=60, res=200):
    """ calculate the periodogram at the specified frequencies, return
    periods, pgram """

    periods = np.linspace(period_low, period_high, res)
    # periods = np.logspace(np.log10(period_low), np.log10(period_high),
    #                       res)
    freqs = 2*np.pi/periods
    pgram = signal.lombscargle(x, y, freqs, precenter=True)

    # significance (see press 1994 numerical recipes, p576)
    var = np.var(y)
    pgram_norm_press = pgram/var
    significance =  1-(1-np.exp(-pgram_norm_press))**len(x)

    #take the normalized power
    pgram_norm = pgram *2 / (len(x) * var)
    return periods, pgram_norm, significance

      