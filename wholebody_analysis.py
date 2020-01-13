# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 13:15:22 2019

@authors: AbelJ ShanY

This code is a basic analysis of whole-body bioluminescence data.
"""
from __future__ import division

import numpy as np
import scipy as sp
import pandas as pd
import dateutil
import minepy
from datetime import datetime
from scipy import interpolate, signal, stats

from matplotlib import gridspec
import matplotlib.pyplot as plt
import prettyplotlib as ppl
from prettyplotlib import brewer2mpl

from LocalImports import PlotOptions as plo
from LocalImports import DecayingSinusoid as ds
from LocalImports import WholeBodyRecording as wbr


# setup the data
INPUT_DIR = 'Demo/WholeBody/M1_030119P2R2GALBWT_Piper1/'

# set up the temperature and humidity data
IMAGING_START = '2019-03-12 17:08:01' #time stamp in format like 19/03/01 17:19:30
IMAGING_INTERVAL = '2 min'
G_FILE = '031219_pump_light_green_liver.xls' # green
R_FILE = '031219_pump_light_red_body.xls' # red

# and the temp and humidity file
TH_FILE = '190301_P2R2GALBWT_Temperature_humidity.xls'
# and actogram file
ACTOGRAM_FILE = '190301_ALB_P2R2G_homo_water_Actogram_Graph Data.csv'

# get the mouse set up
mouse1 = wbr.WholeBodyRecording(INPUT_DIR, R_FILE, G_FILE, IMAGING_START, IMAGING_INTERVAL)

# cut the light intervals
intervals_1 = []
for day in [12,13,14,15,16,17,18]:
    dd = str(day)
    intervals_1 = intervals_1 + [# 12 pulses per day
                 ['2019-03-'+dd+' 07:52:00', '2019-03-'+dd+' 08:24:00'],
                 ['2019-03-'+dd+' 08:52:00', '2019-03-'+dd+' 09:24:00'],
                 ['2019-03-'+dd+' 09:52:00', '2019-03-'+dd+' 10:24:00'],
                 ['2019-03-'+dd+' 10:52:00', '2019-03-'+dd+' 11:24:00'],
                 ['2019-03-'+dd+' 11:52:00', '2019-03-'+dd+' 12:24:00'],
                 ['2019-03-'+dd+' 12:52:00', '2019-03-'+dd+' 13:24:00'],
                 ['2019-03-'+dd+' 13:52:00', '2019-03-'+dd+' 14:24:00'],
                 ['2019-03-'+dd+' 14:52:00', '2019-03-'+dd+' 15:24:00'],
                 ['2019-03-'+dd+' 15:52:00', '2019-03-'+dd+' 16:24:00'],
                 ['2019-03-'+dd+' 16:52:00', '2019-03-'+dd+' 17:24:00'],
                 ['2019-03-'+dd+' 17:52:00', '2019-03-'+dd+' 18:24:00'],
                 ['2019-03-'+dd+' 18:52:00', '2019-03-'+dd+' 19:24:00']
                 ]

for day in [19,20,21,22,23,24]:
    dd = str(day)
    intervals_1 = intervals_1 +[# 4 pulses per day
                 ['2019-03-'+dd+' 07:44:00', '2019-03-'+dd+' 08:32:00'],
                 ['2019-03-'+dd+' 10:44:00', '2019-03-'+dd+' 11:32:00'],
                 ['2019-03-'+dd+' 13:44:00', '2019-03-'+dd+' 14:32:00'],
                 ['2019-03-'+dd+' 16:44:00', '2019-03-'+dd+' 17:32:00']
                 ]

# remove imaging times where lights were on
mouse1.excise_imaging_times(intervals_1, cooldown_ext=5)

# add temp/humidity file
mouse1.import_temperature_humidity(TH_FILE)

# add the activity file
mouse1.import_actogram(ACTOGRAM_FILE)

# perform the processing of imaging data - use the restricted values (i.e., those without the times where lighting was active)
mouse1.process_imaging_data('t_excised','red_excised','green_excised')

# perform the processing of temp/hum data
mouse1.process_temp_hum_data()

# perform the processing of activity data, with a 15 min bin
mouse1.process_activity_data(binsize=15)

# perform the cwt for getting phase info
mouse1.continuous_wavelet_transform()

# perform the correlational analysis for times 24-120, and 216-end
mouse1.correlate_signals(['activity', 'biolum', 'TH'], ['_b', '_es', '_es'], 
                         tmin=110, tmax=182, name='start')
mouse1.correlate_signals(['activity', 'biolum', 'TH'], ['_b', '_es', '_es'], 
                         tmin=216, name='end')


# plot of phase change over time
plt.figure()
plt.plot(mouse1.cwt['activity_es']['x'], np.cos(mouse1.cwt['activity_es']['phase']), 'k', label='Activity')
plt.plot(mouse1.cwt['red_excised_es']['x'], np.cos(mouse1.cwt['red_excised_es']['phase']), 'r:', label='Red - Muscle')
# only plot where green is circadian
good_green = np.where(np.logical_and(mouse1.cwt['green_excised_es']['period']>20, 
                      mouse1.cwt['green_excised_es']['period']<28))[0]
plt.plot(mouse1.cwt['green_excised_es']['x'][good_green], np.cos(mouse1.cwt['green_excised_es']['phase'][good_green]), 'g--', label='Green - Liver')
plt.xlabel('Time (h)')
plt.xticks([0,24,48,72,96,120,144,168,192,216,240,264,288])
plt.ylabel(r'Cos($\phi$), Continuous Wavelet Transform')
plt.ylim([-1,1.5])
plt.legend()
plt.tight_layout(**plo.layout_pad)


# and finally, examine different signals are correlated at the start and at the end of the experiment
fig = plt.figure()
bx = plt.subplot()
pvals = mouse1.corrmat['start']['ps']
rs = mouse1.corrmat['start']['corr']
names = mouse1.corrmat['start']['names']
red_purple = brewer2mpl.get_map('RdBu', 'Diverging', 9).mpl_colormap
ppl.pcolormesh(fig, bx, rs, cmap=red_purple, 
               xticklabels = names, yticklabels = names, vmax=1, vmin=-1)
bx.xaxis.tick_top()
plt.xticks(rotation=90)
# for i in range(len(names)):
#     for j in range(len(names)):
#         if pvals[i, j] < 0.01:
#             text = bx.text(j+0.1, i+1., '**', color="k")
#         elif pvals[i, j] < 0.05:
#             text = bx.text(j+0.25, i+1., '*', color="k")
for tic in bx.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
plt.gca().invert_yaxis()
plt.tight_layout(**plo.layout_pad)


fig = plt.figure()
bx = plt.subplot()
pvals = mouse1.corrmat['end']['ps']
rs = mouse1.corrmat['end']['corr']
names = mouse1.corrmat['end']['names']
red_purple = brewer2mpl.get_map('RdBu', 'Diverging', 9).mpl_colormap
ppl.pcolormesh(fig, bx, rs, cmap=red_purple, 
               xticklabels = names, yticklabels = names, vmax=1, vmin=-1)
# for i in range(len(names)):
#     for j in range(len(names)):
#         if pvals[i, j] < 0.01:
#             text = bx.text(j+0.1, i+1., '**', color="k")
#         elif pvals[i, j] < 0.05:
#             text = bx.text(j+0.25, i+1., '*', color="k")
bx.xaxis.tick_top()
plt.xticks(rotation=90)
for tic in bx.xaxis.get_major_ticks():
    tic.tick1On = tic.tick2On = False
plt.gca().invert_yaxis()
plt.tight_layout(**plo.layout_pad)