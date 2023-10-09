#!/usr/bin/env python
# coding: utf-8

# Seperates the date into on and off files

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

# Add argumenent parsers

parser = argparse.ArgumentParser(
        description='Produces lightcurves and spectra for a set of input FITS files'
        )

parser.add_argument('cropdir', help="The path where the cropped FITS files are. All FITS files (any file ending in .fits) in this directory will be used.")
parser.add_argument('compdir', help="The path where the component FITS files are. Only FITS files that match the file names in CROPDIR will be used, where \"matching\" means the same file name with \"_comp\" inserted before the final \".fits\".")
parser.add_argument('ra', type=float, help="The RA of the source in degrees")
parser.add_argument('dec', type=float, help="The Dec of the source in degrees")
parser.add_argument('--radius', type=float, default=40, help="Only include objects within this angular distance (in arcseconds) from the expected source position (RA, Dec) (default = 30 arcseconds)")
parser.add_argument('--lightcurve_plot', help="The name of the file to which the spectrum plot is saved")
parser.add_argument('--source_name', help="Adding the name of the source to the plots")
parser.add_argument('--alpha', type=float, help="Insert the spectral index")
args = parser.parse_args()

#-----------------------------------------
crop_files = sorted(glob.glob(args.cropdir + "/*.fits"))
basenames = [os.path.basename(f) for f in crop_files]
comp_files = [args.compdir + "/" + f[:-5] + "_comp.fits" for f in basenames]

c = SkyCoord(ra=args.ra*u.degree, dec=args.dec*u.degree)

# Create an empty list to store dictionaries
source_list = []

# Loop over the basenames and corresponding files
for i in range(len(basenames)):
    crop_file = crop_files[i]
    comp_file = comp_files[i]

    # Open FITS files and extract data
    hdul = fits.open(crop_file)
    hdr = hdul[0].header
    date = hdr['DATE-OBS']
    freq = hdr['CRVAL3']
    
    try:
        survey = hdr['PROJECT']
    except KeyError:
        try:
            survey = hdr['SURVEY']
        except KeyError:
            survey = 'unknown'

    df = fits.open(comp_file)
    data = df[1].data
    ra = df[1].data['ra']
    dec = df[1].data['dec']
    fluxes = df[1].data['peak_flux']
    flux_errs = df[1].data['err_peak_flux']
    noise = df[1].data['local_rms']

    if len(ra) != len(fluxes):
        continue     

    # Create coordinates
    catalog = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    # Calculate angular separations    
    seps = c.separation(catalog).to(u.arcsec)

    # Select objects within some user defined (angular) radius
    idx = np.where(seps < args.radius*u.arcsec)[0]

    if len(idx) > 0:
        # Loop over the matching sources and create a dictionary for each
        for j in idx:
            source_dict = {
                'file':basenames[i],
                'date': date,
                'flux': fluxes[j],
                'flux_err': flux_errs[j],
                'freq': freq,
                'noise': noise[j],
                'survey': survey  
            }
            # Append the dictionary to the list
            source_list.append(source_dict)

# Convert the list of dictionaries to a Pandas DataFrame
df = pd.DataFrame(source_list)
print(len(df))
date = np.array(df['date'].values)
print(date)

# Convert DATE to MJD
t = Time(date.astype(str), format='isot', scale='utc')
mjd = t.mjd

# add MJD column to dataframe
df['mjd'] = mjd


# Define the colors and name for each survey
colors = {"TGSS ADR1": "red", "AS107": "blue", "AS207":"orange","AS108": "cyan","AS110":"olive","AS101":"brown","AS113":"PINK","NVSS": "green", "SUMSS":"purple","unknown": "gray"}
survey_name = {"TGSS ADR1": "TGSS", "AS107": "VAST Pilot", "AS207": "VAST","AS108": "GASKAP-HI Pilot", "AS110": "RACS","AS113":"Other ASKAP","AS101":"EMU Pilot","NVSS": "NVSS", "SUMSS":"SUMSS","unknown": "Unknown"}

# Create a new column with the survey names
df["survey_name"] = df["survey"].map(survey_name)

# Create a new column with the survey colors
df["color"] = df["survey"].map(colors)

# Scale the reference frequency 
df['flux_ref']= df['flux']*((887491000/df['freq'])**args.alpha)
df['flux_err_ref']= df['flux_err']*((887491000/df['freq'])**args.alpha)

# Create a new column with SNR values
df['SNR'] = df['flux'] / df['flux_err']

df.to_csv(args.source_name + '_all.csv', index=False)

print(df)

# create empty data frames
on = pd.DataFrame()
off = pd.DataFrame()

# append data to on or off data frames based on SNR value
on = df[df['SNR'] > 5]
off = df[df['SNR'] <= 5]

on.to_csv(args.source_name + '_887000000_on.csv', index=False)
off.to_csv(args.source_name + '_887000000_off.csv', index=False)


