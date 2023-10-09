#!/usr/bin/env python
# coding: utf-8

# Plotting lightcurve and spectrum 

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
from collections import defaultdict
import sys
from scipy.optimize import leastsq
from matplotlib.ticker import FormatStrFormatter
import csv

# Add argumenent parsers

parser = argparse.ArgumentParser(
        description='Produces lightcurves and spectra for a set of input FITS files'
        )

parser.add_argument('cropdir', help="The path where the cropped FITS files are. All FITS files (any file ending in .fits) in this directory will be used.")
parser.add_argument('compdir', help="The path where the component FITS files are. Only FITS files that match the file names in CROPDIR will be used, where \"matching\" means the same file name with \"_comp\" inserted before the final \".fits\".")
parser.add_argument('ra', type=float, help="The RA of the source in degrees")
parser.add_argument('dec', type=float, help="The Dec of the source in degrees")
parser.add_argument('--radius', type=float, default=40, help="Only include objects within this angular distance (in arcseconds) from the expected source position (RA, Dec) (default = 30 arcseconds)")
parser.add_argument('--spectrum_plot', help="The name of the file to which the spectrum plot is saved")
parser.add_argument('--lightcurve_plot', help="The name of the file to which the spectrum plot is saved")
parser.add_argument('--source_name', help="Adding the name of the source to the plots")
args = parser.parse_args()

#-----------------------------------------

crop_files = sorted(glob.glob(args.cropdir + "/*.fits"))
basenames = [os.path.basename(f) for f in crop_files]
comp_files = [args.compdir + "/" + f[:-5] + "_comp.fits" for f in basenames]

c = SkyCoord(ra=args.ra*u.degree, dec=args.dec*u.degree)

output = []
for i in range(len(basenames)):
    crop_file = crop_files[i]
    comp_file = comp_files[i]
    hdul=fits.open(crop_file)
    hdr=hdul[0].header
    dates=hdr['DATE-OBS']
    freq=hdr['CRVAL3']
    
    print('number files:', len(crop_files))

    try:
        survey = hdr['PROJECT']
    except KeyError:
        try:
            survey = hdr['SURVEY']
        except KeyError:
            survey = 'unknown'
    
    df=fits.open(comp_file)
    data=df[1].data
    ra=df[1].data['ra']
    dec=df[1].data['dec']
    fluxes=df[1].data['peak_flux']
    flux_errs=df[1].data['err_peak_flux']
   
    
    if len(ra) != len(fluxes):
        continue

    # Create coordinates
    catalog = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)

    # Calculate angular separations
    seps = c.separation(catalog).to(u.arcsec)

    # Select objects within some user defined (angular) radius
    idx = np.where(seps < args.radius*u.arcsec)[0]


    if len(idx) > 0:
        dates = [dates]*len(idx)
        freq=[freq]*len(idx)
        survey=[survey]*len(idx)
        output += list(zip(dates, fluxes[idx], flux_errs[idx],freq,survey,[basenames[i]]*len(idx)))


#Extract the data from output
date_obs=[row[0] for row in output]
flux=[row[1] for row in output]
flux_err=[row[2] for row in output]
freq=[row[3] for row in output]
survey=[row[4] for row in output]
basenames=[row[5] for row in output]
print('Number Flux:',len(flux))
#Convert time to a Modified Julian Date
t = Time(date_obs, format="isot", scale="utc")
m = t.mjd

# Set the threshold value
SNR = [flux[i] / flux_err[i] for i in range(len(flux))]
on_date = [date_obs[i] for i in range(len(SNR)) if SNR[i] > 5]
on_flux = [flux[i] for i in range(len(SNR)) if SNR[i] > 5]
on_freq = [freq[i] for i in range(len(SNR)) if SNR[i] > 5]
on_flux_err = [flux_err[i] for i in range(len(SNR)) if SNR[i] > 5]
on_survey= [survey[i] for i in range(len(SNR)) if SNR[i] > 5]
on_basenames= [basenames[i] for i in range (len(SNR)) if SNR[i] >5]

output1=list(zip(on_date,on_flux,on_freq,on_flux_err,on_survey))

on_freq = np.array(on_freq)
on_flux = np.array(on_flux)
on_flux_err=np.array(on_flux_err)
print('on_freq', on_freq)
filename = args.source_name + "_5sigma.csv"
header = ["File Name", "Survey", "Date Observed", "On Frequency", "On Flux", "On Flux Error"]
save_data= zip(on_basenames, on_survey, on_date, on_freq, on_flux, on_flux_err)

with open(filename, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(header)
    writer.writerows(save_data)

# spectrum fitting-code
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.size": 7,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

#    "axes.unicode_minus" : False,
#    "pdf.fonttype" : 42, 
#matplotlib.rcParams['ps.fonttype'] = 42
cm = 1/2.54  # centimeters in inches
# http://scipy-cookbook.readthedocs.org/items/FittingData.html
# Define function for calculating a power law
powerlaw = lambda x, amp, index: amp * (x**index)
# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err


def fit_spectrum(on_freq,on_flux,on_flux_err):
    pinit = [10.0, -0.7]
    fit = leastsq(errfunc, pinit, args=(on_freq, on_flux, on_flux), full_output=1)
    covar = fit[1]
    P = fit[0]
    residual = errfunc(P,on_freq, on_flux, on_flux_err)
    chi2red = sum(np.power(residual,2))/(len(on_freq)-len(pinit))
    alpha=P[1]
    amp = np.exp(P[0])
# Errors
    if covar is not None:
        err_alpha = np.sqrt(covar[1][1])
        err_amp = np.sqrt(covar[0][0])
    else:
        err_alpha = None
        err_amp = None
    return alpha, err_alpha, amp, err_amp, chi2red

# define the colors for each survey

# define the colors for each survey
colors = {"TGSS ADR1": "red", "AS107": "blue", "AS207":"orange","AS108": "cyan","AS110":"olive","AS101":"brown","AS113":"PINK","NVSS": "green", "SUMSS":"purple","unknown": "gray"}
new_names = {"TGSS ADR1": "TGSS", "AS107": "VAST Pilot", "AS207": "VAST","AS108": "GASKAP-HI Pilot", "AS110": "RACS","AS113":"Other ASKAP","AS101":"EMU Pilot","NVSS": "NVSS", "SUMSS":"SUMSS","unknown": "Unknown"}


# collect the flux, frequency, and flux_err data for each survey
data = defaultdict(lambda: {'on_freq': [], 'on_flux': [], 'on_flux_err': [], 'on_survey': []})

print(on_freq)
print(on_flux)

freqcent = np.mean(on_freq)
alpha, err_alpha, amp, err_amp, chi2red = fit_spectrum(np.log(on_freq/freqcent), np.log(on_flux), on_flux_err/on_flux)

all_freqs = []
all_flux = []
all_flux_err=[]

# Plot the flux vs frequency for each survey
fig, ax = plt.subplots()

for item in output1:
    on_date,on_flux,on_freq,on_flux_err,on_survey = item
    data[on_survey]['on_freq'].append(on_freq)
    data[on_survey]['on_flux'].append(on_flux)
    data[on_survey]['on_flux_err'].append(on_flux_err)
    data[on_survey]['on_survey'].append(on_survey)

#Plot the flux vs frequency for each survey
for on_survey, survey_data in data.items():
    on_freq = survey_data['on_freq']
    on_flux = survey_data['on_flux']
    on_flux_err = survey_data['on_flux_err']
    
    # append all the on_freq values to the list
    all_freqs.extend(on_freq)
    all_flux.extend(on_flux)
    all_flux_err.extend(on_flux_err)

    #Check if survey is in the colors dictionary and set color and label accordingly
    if on_survey in colors:
        color = colors[on_survey]
        label = new_names[on_survey]
    else:
        color = "gray"
        label = on_survey

    on_freq = np.array(on_freq)
    ax.errorbar(on_freq/1.e6, on_flux, yerr= on_flux_err, lw=0.5, elinewidth=0.5, fmt=".", markeredgewidth=0.5, label=label)


all_freq = np.concatenate([survey_data['on_freq'] for survey_data in data.values()])
all_freq = np.unique(all_freq)
all_flux = np.concatenate([survey_data['on_flux'] for survey_data in data.values()])
all_flux = np.unique(all_flux)
all_flux_err = np.concatenate([survey_data['on_flux_err'] for survey_data in data.values()])
all_flux_err = np.unique(all_flux_err)
ax.plot(all_freq/1.e6, powerlaw(all_freq/freqcent, amp, alpha), color="blue", zorder=10, lw=0.5)
#Set the labels and title
ax.set_xlabel('Frequency [MHz]')
ax.set_ylabel('Flux [Jy]')
ax.set_title(args.source_name +' Spectrum')

if err_alpha is None:
    ax.set_title(args.source_name +' Spectrum'+'   '   r"$\alpha = {0:4.2f}$, $\chi^2_r = {1:4.2f}$".format(alpha, chi2red))
else:
    ax.set_title(args.source_name +' Spectrum'+'   '   r"$\alpha = {0:4.2f}\pm{1:4.2f}$, $\chi^2_r = {2:4.2f}$".format(alpha, err_alpha, chi2red))


#Plot it on a log scale 
ax.set_yscale('log')
ax.set_xscale('log')
#Set the tick values and labels
ax.set_yticks([0.001, 0.01,0.1], ['0.001','0.01','0.1'])
ax.set_xticks([100,200,300,400,500,600,1000])
ax.set_xticklabels(['100','200','300','400','500','600', '1000'])
ax.legend()
print(alpha)

#Saving the plot
if args.spectrum_plot is None:
    plt.show()
else:
    fig.savefig(args.spectrum_plot)

	
