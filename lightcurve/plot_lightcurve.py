#!/usr/bin/env python
# coding: utf-8

# Plots the lightcurve

#Import relevant functions
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import argparse
import os
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from collections import defaultdict
import pandas as pd
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec


# Add argumenent parsers

parser = argparse.ArgumentParser(
        description='Produces a lightcurve for the variable source')

parser.add_argument('on', help="The path where the csv file is for the on data for the variable source")
parser.add_argument('off', help="The path where the csv file is for the off data for the refernce source")

parser.add_argument('--lightcurve_plot', help="The name of the file to which the spectrum plot is saved")
parser.add_argument('--source_name', help="Adding the name of the source to the plots")
args = parser.parse_args()

on = pd.read_csv(args.on)
off = pd.read_csv(args.off)

fig, ax = plt.subplots(figsize=(8, 6))

# creating a subplot with no data to display the legend
gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

# plot the varaible source when its off
for survey in on['survey_name'].unique():
        survey_data = on[on['survey_name'] == survey]
        ax1.errorbar(survey_data['mjd'].values, survey_data['flux_ref'].values, yerr=survey_data['flux_err_ref'].values, fmt='.', color=survey_data['color'].iloc[0], label=survey, zorder=2)

# plot the variable source when its off
for survey in off['survey_name'].unique():
        survey_data = off[off['survey_name'] == survey]
        ax1.errorbar(survey_data['mjd'].values, survey_data['flux_ref'].values, yerr=survey_data['flux_err_ref'].values, marker="v", color=survey_data['color'].iloc[0], label=survey, zorder=2, linestyle='none')


# set the labels and title
ax1.set_xlabel('Frequency [Hz]')
ax1.set_xlabel('Date [MJD]')
ax1.set_ylabel('Flux [Jy]')
ax1.set_title(args.source_name + ' lightcurve')

# add legend to ax2
handles, labels = ax1.get_legend_handles_labels()
ax2.legend(handles, labels, loc='center', frameon=False)
ax2.axis('off')

#Saving the plot
if args.lightcurve_plot is None:
    plt.show()
else:
    fig.savefig(args.lightcurve_plot)
