# Signal Processing
After raw photometry data is loaded in to MATLAB, signal processing is conducted to account for photobleaching, motion artifacts, and other sources of "noise". The signal processing functions are written to be as flexible to differing streams and naming conventions as possible, but if the functions don't match your data, please reach out and let us know and we will update.

## Background Scaling and Subtraction
To correct for photobleaching and motion artifact, the background is scaled to the signal stream and subtracted. While previous tools and packages have used regression based approaches, PASTa scales the background based on the frequency domain of the streams, preventing over or underfitting and maintaining the shape of the background. 

### Background Scaling Factor
To determine the scaling factor for the background stream relative to the signal stream, both the signal and background streams are converted to the frequency domain via Fast Fourier Transform (FFT). The scaling factor is calculated as for all frequencies greater than 10 Hz. The background stream is scaled in the time domain; the background is centered around zero, multiplied by the scaling factor, and then adjusted to center around the signal mean.

### Subtraction
The scaled background is subtracted from the raw signal. Subtracted signal is output as dF/F. Users have the ability to override PASTa protocol defaults to modify parameters including the scaling threshold frequency cutoff, subtracted data output, or select an alternative method of background scaling.


### Other Scaling Methods
In our hands, frequency scaling performs better than other scaling approaches. Advantages are particularly notable for sensors such as GRABDA2H, for which 405nm is not a perfect isosbestic control. Use of frequency scaling rescues the use of the 405nm stream to control for photobleaching and motion artifact. This may be
particularly useful as new sensors are continually in development, not all of which have an isosbestic or commercially available control wavelength for use in photometry systems.

However, users may want to compare to other commonly used methods of scaling the background to the signal. The subtractFPdata function includes options to use OLS regression, detrending and OLS regression, Lowess Smoothing and OLS Regression, and IRLS regression to scale the background prior to subtraction.


## Filtering
After subtraction, the signal is filtered to remove high frequency noise and facilitate further analysis. PASTa protocol has default filter settings, but users can override these as required by sensor or experimental design.

* __Filter type:__ defaults to a 3rd order Butterworth band-pass filter, which preserves the shape of the signal in the pass band.

* __Filter cutoffs:__ The high-pass cutoff is set to 0.0051 Hz to remove the zero frequency component of the power spectrum. The low-pass cutoff is set to 2.2860 Hz to preserve the frequencies of interest while attenuating high frequency noise.

Users can override defaults to modify the filter type (band-pass, high-pass, low-pass), order, and cutoff frequencies.

## Normalization
Normalization converts the filtered signal to Z score. Multiple methods are included in the tool box to accomodate a variety of experimental designs.
Whole session uses the mean and SD of the whole session. Session baseline uses the mean and SD from a specified session baseline, which may be useful in cases where a drug is delivered mid session. Data can also be normalized on an individual trial basis to a local pre-trial baseline.

Differing normalization methods may be advantageous depending on the experimental design. If additional options are required, please let us know!

