#!/usr/bin/env python
# coding: utf-8

# Plot lightcurve for the refernce and variable source m

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
        description='Produces a lightcurve plotting both the variable source and reference')

parser.add_argument('on', help="The path where the csv file is for the on data for the variable source")
parser.add_argument('off', help="The path where the csv file is for the off data for the refernce source")

parser.add_argument('on_ref', help="The path where the csv file is for the on data for the reference source")
parser.add_argument('off_ref', help="The path where the csv file is for the off data for the refernce source")

parser.add_argument('--both_lightcurves', help="The name of the file to which the lightcurveplot is saved")
parser.add_argument('--source_name', help="Adding the name of the source to the plots")
args = parser.parse_args()

on = pd.read_csv(args.on)
off = pd.read_csv(args.off)

on_ref = pd.read_csv(args.on_ref)
off_ref = pd.read_csv(args.off_ref)

#Normalise the reference source by their mean

on_ref['flux_ref_norm']=on_ref['flux_ref']/np.mean(on_ref['flux_ref'])
off_ref['flux_ref_norm']=off_ref['flux_ref']/np.mean(off_ref['flux_ref'])

on_ref['flux_err_ref_norm']=on_ref['flux_err_ref']/np.mean(on_ref['flux_ref'])
off_ref['flux_err_ref_norm']=off_ref['flux_err_ref']/np.mean(off_ref['flux_ref'])



# Normalise the reference source and then multiple by 1.5 to seperate from the refernce source

# Calculate the mean flux value
#mean_flux_ref_on = on_ref['flux_ref'].mean()
#mean_flux_ref_off = off_ref['flux_ref'].mean()

# Normalize the flux and error by the mean flux value
#on_ref['flux_ref_norm'] = (on_ref['flux_ref'] / mean_flux_ref_on)*0.03
#on_ref['flux_err_ref_norm'] = (on_ref['flux_err_ref'] / mean_flux_ref_on)*0.03

#off_ref['flux_ref_norm'] = off_ref['flux_ref'] / mean_flux_ref_off
#off_ref['flux_err_ref_norm'] = off_ref['flux_err_ref'] / mean_flux_ref_off


# Change the survey name for the reference source
#survey_name_ref = {"TGSS" : "TGSS Ref", "VAST Pilot" : "VAST Pilot Ref", "VAST" : "VAST Ref","GASKAP-HI Pilot": "GASKAP-HI Pilot Ref", "RACS": "RACS Ref","Other ASKAP": "Other ASKAP Ref", "EMU Pilot": "EMU Pilot Ref","NVSS": "NVSS Ref", "SUMSS":"SUMSS Ref","unknown": "Unknown Ref"}
#color_ref = {"red" : "blue",  "blue" : "green", "orange" : "red", "cyan" : "olive","olive" : "maroon","brown": "black","PINK" : "blue", "green" : "cyan", "purple": "black", "gray":"grey"}

# Use map function to replaxes the survye names with the new names 
#on_ref["survey_name_ref"] = on_ref["survey_name"].map(survey_name_ref)
#off_ref["survey_name_ref"] = off_ref["survey_name"].map(survey_name_ref)

# Use map function to replaxes the survye names with new colors
#on_ref["color_ref"] = on_ref["color"].map(color_ref)
#off_ref["color_ref"] = off_ref["color"].map(color_ref)

#print(on_ref)


fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 2, width_ratios=[5, 1])

ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[1, 0])
ax3 = plt.subplot(gs[1, 1])
ax4 = plt.subplot(gs[0, 1])

print(len(on['flux_ref']))
print(len(off['flux_ref']))
print(len(on_ref['flux_ref']))
print(len(off_ref['flux_ref']))


# plot the varaible source when its off
for survey in on['survey_name'].unique():
	survey_data = on[on['survey_name'] == survey]	
	ax1.errorbar(survey_data['mjd'].values, survey_data['flux_ref'].values, yerr=survey_data['flux_err_ref'].values, fmt='.', color=survey_data['color'].iloc[0], label=survey, zorder=2)

# plot the variable source when its off
for survey in off['survey_name'].unique():
	survey_data = off[off['survey_name'] == survey]
	ax1.errorbar(survey_data['mjd'].values, survey_data['flux_ref'].values, yerr=survey_data['flux_err_ref'].values, marker="v", color=survey_data['color'].iloc[0], label=survey, zorder=2, linestyle='none')

# plot the on survey for the reference source 
for survey in on_ref['survey_name'].unique():
        survey_data = on_ref[on_ref['survey_name'] == survey]
        ax2.errorbar(survey_data['mjd'].values, survey_data['flux_ref'].values, yerr=survey_data['flux_err_ref'].values, fmt='.', color=survey_data['color'].iloc[0], label=survey, zorder=2)

# plot the off survey for the reference source
for survey in off_ref['survey_name'].unique():
        survey_data = off_ref[off_ref['survey_name'] == survey]
        ax2.errorbar(survey_data['mjd'].values, survey_data['flux_ref'].values, yerr=survey_data['flux_err_ref'].values, marker="v", color=survey_data['color'].iloc[0], label=survey, zorder=2, linestyle='none')

#Plot the normalised reference flux by its mean

# plot the on survey for the reference source 
#for survey in on_ref['survey_name'].unique():
 #       survey_data = on_ref[on_ref['survey_name'] == survey]
  #      ax2.errorbar(survey_data['mjd'].values, survey_data['flux_ref_norm'].values, yerr=survey_data['flux_err_ref_norm'].values, fmt='x', color=survey_data['color'].iloc[0], label='{}: Normalised'.format(survey), zorder=2)

# plot the off survey for the reference source
#for survey in off_ref['survey_name'].unique():
 #       survey_data = off_ref[off_ref['survey_name'] == survey]
  #      ax2.errorbar(survey_data['mjd'].values, survey_data['flux_ref_norm'].values, yerr=survey_data['flux_err_ref_norm'].values, marker="x", color=survey_data['color'].iloc[0], label='{}: Normalised'.format(survey), zorder=2, linestyle='none')


# set the labels and title
ax1.set_xlabel('Frequency [Hz]')
ax1.set_xlabel('Date [MJD]')
ax1.set_ylabel('Flux [Jy]')
ax1.text(np.max(on['mjd']), np.max(off['flux_ref'])+0.0005, 'Frequency = 887491000 MHz', ha='right', va='top')

ax2.set_xlabel('Frequency [Hz]')
ax2.set_xlabel('Date [MJD]')
ax2.set_ylabel('Reference Flux [Jy]')



# add legend to ax4
handles, labels = ax1.get_legend_handles_labels()
ax4.legend(handles, labels, loc='center', frameon=False)
ax4.axis('off')

# add legend to ax3
handles, labels = ax2.get_legend_handles_labels()
ax3.legend(handles, labels, loc='center', frameon=False)
ax3.axis('off')


# Add a title for the entire grid space
fig.suptitle(args.source_name + ' lightcurve', fontsize=16)

#Saving the plot
if args.both_lightcurves is None:
    plt.show()
else:
    fig.savefig(args.both_lightcurves)
