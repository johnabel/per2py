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
from LocalImports import ProcessingFunctions as pf

#inputs nms pre
pull_from_imagej = False
input_folder = 'Data/scn2_VIPBmal1KO_20170909_SGLE1/' # edit this
input_ij_extension = '.csv'# edit this

# do we want plots?
plotcount = 2
# for each dataset, this is formatted [descriptor, root filename, color ('Red', or 'Green')]
input_ij_file_type   = [['Pre_Green', '090917VIPMOP_Pre_Green'],
                        ['Pre_Green', '090917VIPMOP_Pre_Red'],
                        ['TTX_Red', '090917VIPMOP_TTX_Red'],
                        ['TTX_Green', '090917VIPMOP_TTX_Green'],
                        ['Wash Red', '090917VIPMOP_Wash_Red'],
                        ['Wash_Green', '090917VIPMOP_Wash_Green']
                        ] # edit this

# list all the datasets
all_inputs=[]
for input_ij in input_ij_file_type:
    all_inputs.append(pf.generate_filenames_dict(input_folder, input_ij[1], 
                                    pull_from_imagej, input_ij[0], 
                                    input_ij_extension, input_ij[2]))