# Lightcurve plotting code for the sources

file_lightcurve.py contains the code to sperate the data into on and off frames

```
python /path/to/file_lightcurve.py /path/to/crop /path/to/comp ra dec --source_name  --alpha
```

plot_lightcurve.py contains the code to plot the lightcurves

```
python /path/to/plot_lightcurve.py /path/to/on /path/to/off --lightcurve_plot --source_name
```

plot_both_lightcurve.py contains the code to plot the lightcurve of the source and the reference source together

```
python /path/to/plot_both_lightcurve.py /path/to/on /path/to/off /path/to/ref_on /path/to/ref_off --source_name  --both_lightcurves
```
