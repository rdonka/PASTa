# Function Documentation Overview
This page contains additional documentation for each function within PASTa, as well as examples of inputs.

# Data Preparation Functions
This set of functions is used to prepare raw photometry data, match it with experimental metadata, and load data into a structure in MATLAB. Functions are provided to handle data collected via TDT equipment and software Synapse, or a generic file structure with data streams saved to CSV files.

## loadKeys
Combines subject key and file key into a data structure, and appends the provided computeruserpath to the paths in the file key.

**INPUTS:**

* __COMPUTERUSERPATH:__ A variable containing the unique portion of the file explorer path for the users specific computer. For example, 'C:\Users\rmdon\'. Make sure the computeruserpath ends in a forward slash.
* __SUBJECTKEYNAME:__ A variable containing a string with the name of the subject key csv file for the experiment (see _"3. Data Analysis Pipeline"_ for more details). To omit a subject key and only load in the file key, set subjectkeyname = "".
* __FILEKEYNAME:__ A variable containing a string with the name of the file key csv file for the experiment (see _"3. Data Analysis Pipeline"_ for more details).

**OUTPUTS:**

* __EXPERIMENTKEY:__ A data structure called "experimentkey" that includes the joined file key and subject key with the computer user path appended to raw and extracted folder paths.

**EXAMPLE:**
```
computeruserpath =  'C:\Users\MYNAME\'; % Computer specific portion of file navigation paths
subjectkeyname = 'Subject Key.csv'; % Name of csv file containing subject information; set to '' if not using a subject key
filekeyname = 'File Key.csv'; % Name of csv file containing session information, raw data folder names, and paths

[experimentkey] = loadKeys(computeruserpath, subjectkeyname, filekeyname); % Load keys into a data structure called experimentkey
```

**NOTES:**

* FILEKEY must contain at a minimum the fields _Subject_, _RawFolderPath_, and _ExtractedFolderPath_. 
* SUBJECTKEY must contain at a minimum the field _Subject_
* Folder paths must end with a slash.
* The subject and file keys are joined based on Subject ID. Subject key must contain every subjects in the file key. If there is a mismatch, you will receive an error message that the right table does not contain all the key variables that are in the left table. The error message will display the unique subject IDs present in each key so you can determine where the mismatch occurred.
* Fields in subject and file key must be named uniquely. The only field that should be named the same in both keys is Subject.

## extractTDTdata
This function is used to extract TDT data from saved blocks recorded via the software _Synapse_. For each block, extractTDTdata calls the function "TDTbin2mat" (TDT, 2019) and inputs the RawFolderPath to extract fiber photometry data recorded with Synapse. Extracted blocks are parsed it into a single data structure containing all fields, streams, and epocs. The function will identify the signal channel by matching the names in the input SIGSTREAMNAMES and the control channel by matching the names in the input BAQSTREAMNAMES. The name inputs can include a list of stream names if channel naming conventions vary by rig. Each block is saved as a separate data structure in a '.mat' file at the location specified by the inputs in extractedfolderpaths. 

**INPUTS:**

* __RAWFOLDERPATHS:__ a string array containing the paths to the folder location of the raw data blocks to be extracted. The string array should contain one column with each full path in a separate row.
* __EXTRACTEDFOLDERPATHS:__ a string array containing the paths to the folder location in which to save the extracted MatLab structs for each block to be extracted. The string array should contain one column with each full path in a separate row.  
* __SIGSTREAMNAMES:__ A cell array containing strings with the names of the streams to be treated as signal. Note that only one stream per file can be treated as signal. If different files have different stream names, include all stream names in the cell array.
* __BAQSTREAMNAMES:__ A cell array containing strings with the names of the streams to be treated as background. Note that only one stream per file can be treated as background. If different files have different stream names, include all stream names in the cell array.

**OPTIONAL INPUTS:**

* __CLIP:__ the number of seconds to remove on either end of the data streams. If not specified, defaults to 5 seconds.
* __SKIPEXISTING:__ A binary variable containing a 0 if pre-existing extracted blocks should be re-extracted or a 1 if pre-existing extracted blocks should be skipped. This allows the user to toggle whether or not to extract every block, or only blocks that have not previously been extracted. If not specified, defaults to 1 (skip previously extracted blocks).

**OUTPUTS:**

Saved .mat data structures for each block in the location specified by extractedfolderpaths.

**EXAMPLE - DEFAULT:**
```
sigstreamnames = {'x65A', '465A'}; % All names of signal streams across files
baqstreamnames = {'x05A', '405A'}; % All names of background streams across files
rawfolderpaths = string({experimentkey.RawFolderPath})'; % Create string array of raw folder paths
extractedfolderpaths = string({experimentkey.ExtractedFolderPath})'; % Create string array of extracted folder paths

extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames,clip,skipexisting); % extract data
```

**EXAMPLE - MANUALLY SPECIFIED CLIP AND SKIPEXISTING:**
```
clip = 3;
skipexisting = 0;
sigstreamnames = {'x65A', '465A'}; % All names of signal streams across files
baqstreamnames = {'x05A', '405A'}; % All names of background streams across files
rawfolderpaths = string({experimentkey.RawFolderPath})'; % Create string array of raw folder paths'
extractedfolderpaths = string({experimentkey.ExtractedFolderPath})'; % Create string array of extracted folder paths'

extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames,'clip',clip,'skipexisting',skipexisting); % extract data
```

## loadTDTdata
For use with previously extracted data collect with TDT equipment and software Synapse. loadTDTData loads previously extracted .mat data blocks into a data structure for further analysis. Each block is one row. The input experiment key must include at a minimum the extracted folder path location of the block. Any additional information about the subject and session in the experiment key will be matched to the extracted data.

**INPUTS:**

* __EXPERIMENTKEY:__ A prepared data structure with at minimum the ExtractedFolderPath containing the full path to the individual session data structures to be loaded. The EXPERIMENTKEY can be prepared with the LOADKEYS function.

**OUTPUTS:**

* __DATA:__ the input data structure with each individual extracted block added by row.

**EXAMPLE:**
```
[rawdata] = loadKeydata(experimentkey); % Load data based on the experiment key into the structure 'rawdata'
```

## loadCSVdata
For use with data collected and stored to a general file structure. 

_coming soon_


## trimFPdata
Used to trim samples at the very start and end of recordings that are not to be included in analysis (such as the the first two minutes of the session, or the first samples before a hardware control program is initiated). Trims all specified data streams from the index in trimstart to the index in trimend, and adjusts epocs by the amount trimmed by trimstart. Users must pre-prepare the trim start and end indexes to specify as inputs for the function.

**INPUTS:**

* __DATA:__ A data frame containing at least the specified input fields.
* __TRIMSTART:__ The location to start trimming at.
* __TRIMEND:__ The location to end trimming at.
* __WHICHSTREAMS:__ A cell array containing the names of all the streams to be trimmed.

**OPTIONAL INPUTS:**

* __WHICHEPOCS:__ A cell array containing the names of all the epocs to be adjusted due to trimming - subtract the (start loc - 1) from the epoc.

**OUTPUTS:**

* __DATA:__ The data structure with the specified data stream containing the trimmed data.

**EXAMPLE:**
```
trimstart = 'sessionstart'; % name of field with session start index
trimend = 'sessionend'; % name of field with session end index
whichstreams = {'sig', 'baq','time'}; % which streams to trim
whichepocs = {'injt','sess'}; % which epocs to adjust to maintain relative position

[data] = trimFPdata(rawdata,trimstart,trimend, whichstreams,whichepocs); % Output trimmed data into new structure called data
```

# Signal Processing Functions

## subtractFPdata
Used to subtract the background photometry stream (eg, 405nm) from the signal stream (eg, 465nm), convert the subtracted signal to delta F/f, and apply a filter to denoise the output. Users must input the data structure with the raw data, the names of the fields containing the signal and background streams, and the sampling rate of the collected data.

**INPUTS:**

* __DATA:__ A data frame containing at least the specified input fields.
* __SIGFIELD:__ The name (string) of the field containing the signal stream.
* __BAQFIELD:__ The name (string) of the field containing the background stream.
* __WHICHFS:__ The name (string) of the field containing the sampling rate of the raw data collection in Hz.

**OPTIONAL INPUTS:**

* __BAQSCALINGTYPE:__ A string to specify the type of background scaling to apply. Options are 'frequency', 'sigmean', 'OLS', 'detrendOLS', 'smoothedOLS', or 'IRLS'. Default: 'frequency'.
    * _'frequency':_ Scales the background to the signal channel based on ratio of specified frequency bands in the FFT (frequency domain) of the channels.
    * _'sigmean':_ Scales the background to the signal channel based on the ratio of the mean of the signal to the mean of the background (time domain).
    * _'OLS':_ Uses ordinary least-squares regression to generate scaled background.
    * _'detrendOLS':_ Removes the linear trend from signal and background streams prior to using ordinary least-squares regression to generate scaled background.
    * _'smoothedOLS':_ Applies lowess smoothing to the signal and background streams prior to using ordinary least-squares regression to generate scaled background.
    * _'IRLS':_ Uses iteratively reweighted least squares regression to generate scaled background.
* __BAQSCALINGFREQ:__ Only used with 'frequency' scaling. Numeric frequency (Hz) threshold for scaling the background to signal channel. Frequencies above this value will be included in the scaling factor determination. Default: 10 Hz.
* __BAQSCALINGPERC:__ Only used with 'frequency' and 'sigmean' scaling. Adjusts the background scaling factor to be a percent of the derived scaling factor value. Default: 1 (100%).
* __SUBTRACTIONOUTPUT:__ Output type for the subtracted data. Default: 'dff'
    * _'dff':_ Outputs subtracted signal as delta F/F.
    * _'df':_ Outputs subtracted signal as delta F.
* __FILTERTYPE:__ A string to specify the type of filter to apply after subtraction. Default: 'bandpass'.
    * _'nofilter':_ No filter will be applied.
    * _'bandpass':_ A bandpass filter will be applied.
    * _'highpass':_ Only the high pass filter will be applied.
    * _'lowpass':_ Only the low pass filter will be applied.
* __PADDING:__ Defaults to 1, which applies padding. Padding takes the first 10% of the stream, flips it, and appends it to the data before filtering. Appended data is trimmed after filtration. Set to 0 to turn off padding of data streams. Default: 1.
* __PADDINGPERC:__ Percent of data length to use to determine the number of samples to be appended to the beginning and end of data in padding. Set to minimum 10%. Default: 0.1 (10%).
* __FILTERORDER:__  The order to be used for the chosen butterworth filter. Default: 3.
* __HIGHPASSCUTOFF:__  The cutoff frequency (hz) to be used for the high pass butterworth filter. Default: 2.2860.
* __LOWPASSCUTOFF:__ The cutoff to be used for the low pass butterworth filter. Default: 0.0051.
    * NOTE: 'bandpass' applies both the high and low cutoffs to design the filter.
* __SUPRESSDISP:__ If set to anything other than 0, this will suppress the command window displays. Default: 0.

**OUTPUTS:**

* __DATA:__ The original data structure with added fields with the scaled background ('baq_scaled'), subtracted signal ('sigsub'), and subtracted and filtered signal ('sigfilt'). All inputs and defaults will be added to the data structure under the field 'inputs'.
    * __NOTE:__ If using BAQSCALINGMETHOD _'detrendOLS'_, additional fields containing the detrended signal and background ('sig_detrend' and 'baq_detrend') will be added to the data frame. If using BAQSCALINGMETHOD _'smoothedOLS'_, additional fields containing the smoothed signal and background ('sig_smoothed' and 'baq_smoothed') will be added to the data frame.


**EXAMPLE - DEFAULT:**
```
sigfield = 'sig';
baqfield = 'baq';
fs = 1017;

data = subtractFPdata(data,sigfield,baqfield,fs);
```

**EXAMPLE - Frequency Scaling with 20Hz Threshold and Highpass Filter Only:**
```
sigfield = 'sig';
baqfield = 'baq';
fs = 1017;

data = subtractFPdata(data,sigfield,baqfield,fs,'baqscalingfreq',20,'filtertype,'highpass');
```


## normSession
Used to normalize the subtracted and filtered signal stream to Z score based on the whole session mean and standard deviation.

**INPUTS:**

* __DATA:__ A data frame containing at least the stream specified to be normalized.
* __WHICHSTREAM:__ A string with the name of the field containing the stream to be normalized.


**OUTPUTS:**

* __DATA:__ The data structure with the field 'data.WHICHSTREAMz_normsession' containing the whole session normalized signal added.

**EXAMPLE:**
```
[data] = normSession(data,'sigfilt'); % Outputs whole session z score

```

## normBaseline
Used to normalize the subtracted and filtered signal stream to Z score based on the a session baseline mean and standard deviation.

**INPUTS:**

* __DATA:__ A data frame containing at least the stream specified to be normalized.
* __WHICHSTREAM:__ A string with the name of the field containing the stream to be normalized.
* __WHICHBLSTART:__ A string with the name of the field containing the index of the start of the baseline period.
* __WHICHBLEND:__ A string with the name of the field containing the index of the end of the baseline period.


**OUTPUTS:**

* __DATA:__ The data structure with the field 'data.WHICHSTREAMz_normbaseline' containing the full stream normalized to baseline added.

**EXAMPLE:**
```
% Prepare baseline start and end indices
for eachfile = 1:length(data) % prepare indexes for baseline period start and end
    data(eachfile).BLstart = 1;
    data(eachfile).BLend =  data(eachfile).injt(1);
end

% Normalize filtered signal to session baseline mean and SD
[data] = normBaseline(data,'sigfilt','BLstart','BLend');
```


# Transient Detection and Quantification Functions
## findSessionTransients
Main function to detect and quantify transient events in whole session data streams. Use of this main function will set default values for parameters not specific as optional inputs. Users must input the raw data structure containing the data stream, a field with threshold values for transient inclusion, and the sampling rate of the stream. Optional inputs can be specified to modify transient detection parameters, but defaults will be applied if not specified.


**INPUTS:**

* __DATA:__ This is a structure that contains at least the data stream you want to analyze, a field with the threshold values, and a field with the sampling rate of the data stream.
* __WHICHBLTYPE:__ A string with the type of pre-transient baseline to use for amplitude inclusion and quantification.
    * _'blmean':_ Pre-transient baselines are set to the mean of the pre-transient window.
    * _'blmin':_ Pre-transient baselines are set to the minimum value within the pre-transient window.
    * _'localmin':_ Pre-transient baselines are set to the local minimum directly preceding the transient within the baseline window.
* __WHICHSTREAM:__ The name (string) of the field containing the stream to be analyzed for transients.
* __WHICHTHRESHOLD:__ The name (string) of the field containing the prepared numeric threshold values for each stream. For example, 'threshold_3SD'. 
* __WHICHFS:__ The name (string) of the field containing the sampling rate of the raw data collection in Hz.


**OPTIONAL INPUTS:**

* __PREMINSTARTMS:__ Number of millseconds pre-transient to use as the start of the baseline window. _Default: 800_
* __PREMINENDMS:__ Number of millseconds pre-transient to use as the end of the baseline window. _Default: 100_
* __POSTTRANSIENTMS:__ Number of millseconds post-transient to use for the post peak baseline and trimmed data output. _Default: 2000_
* __QUANTIFICATIONHEIGHT:__ The height at which to characterize rise time, fall time, peak width, and AUC. Must be a number between 0 and 1. _Default: 0.5_
* __OUTPUTTRANSIENTDATA:__ Set to 1 to output cut data streams for each transient event for further analysis or plotting. Set to 0 to skip. _Default: 1_

**OUTPUTS:**

* __DATA:__ The original data structure with sessiontransients_WHICHBLTYPE_THRESHOLDLABEL added in. For more details on individual variables, see the PASTa user guide page _Transient Detection_. The output contains four nested tables: 
    * _INPUTS:_ Includes all required and optional inputs. If optional inputs are not specified, defaults will be applied.
    * _TRANSIENTQUANTIFICATION:_ Includes the quantified variables for each transient, including amplitude, rise time, fall time, width, and AUC. 
    * _TRANSIENTSTREAMLOCS:_ Pre-transient baseline, transient peak, rise, and fall locations for each transient to match the cut transient stream data.
    * _TRANSIENTSTREAMDATA: Cut data stream from baseline start to the end of the post-transient period for each transient event.


**NOTES:**

* Threshold values should be calculated before using the findSessionTransients functions. Typically thresholds are set to 2-3 SDs. If the input data stream is Z scored, this can be the actual SD threshold number. If the input data stream is not Z scored, find the corresponding value to 2-3 SDs for each subject.
* For transient data outputs, each transient is in a separate row. 
* If OUTPUTTRANSIENTDATA is set to anything other than 1, the TRANSIENTSTREAMLOCS and TRANSIENTSTREAMDATA tables will be skipped and not included in the output.


**EXAMPLE:**
```
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients based on pre-peak baseline window minimum
[data] = findSessionTransients(data,'blmin','sigfiltz_normsession_trimmed','threshold3SD','fs');

% Find session transients based on pre-peak baseline window mean
[data] = findSessionTransients(data,'blmean','sigfiltz_normsession_trimmed','threshold3SD','fs','preminstartms',600);

% Find session transients based on pre-peak local minimum (last minumum before the peak in the baseline window)
[data] = findSessionTransients(data,'localmin','sigfiltz_normsession_trimmed','threshold3SD','fs','quantificationheight',0.25);
```


### findSessionTransients_blmean
Sub function to detect and quantify transient events in whole session data streams. Transient baselines are defined as the mean of the pre-transient baseline window.

This function is called by the main function _findSessionTransients_, which sets default parameters for transient detection and quantification. If called outside the main function, users must specify all input values manually.

**INPUTS:**

* __DATA:__ This is a structure that contains at least the data stream you want to analyze, a field with the threshold values, and a field with the sampling rate of the data stream.
* __WHICHSTREAM:__ The name (string) of the field containing the stream to be analyzed for transients.
* __WHICHTHRESHOLD:__ The name (string) of the field containing the prepared numeric threshold values for each stream. For example, 'threshold_3SD'. 
* __WHICHFS:__ The name (string) of the field containing the sampling rate of the raw data collection in Hz.
* __PREMINSTARTMS:__ Number of millseconds pre-transient to use as the start of the baseline window. _Reccomended value: 800_
* __PREMINENDMS:__ Number of millseconds pre-transient to use as the end of the baseline window. _Reccomended value: 100_
* __POSTTRANSIENTMS:__ Number of millseconds post-transient to use for the post peak baseline and trimmed data output. _Reccomended value: 2000_
* __QUANTIFICATIONHEIGHT:__ The height at which to characterize rise time, fall time, peak width, and AUC. Must be a number between 0 and 1. _Reccomended value: 0.5_
* __OUTPUTTRANSIENTDATA:__ Set to 1 to output cut data streams for each transient event for further analysis or plotting. Set to 0 to skip. _Reccomended value: 1_

**OUTPUTS:**

* __DATA:__ The original data structure with sessiontransients_blmean_THRESHOLDLABEL added in. For more details on individual variables, see the PASTa user guide page _Transient Detection_. The output contains four nested tables: 
    * _INPUTS:_ Includes all function inputs.
    * _TRANSIENTQUANTIFICATION:_ Includes the quantified variables for each transient, including amplitude, rise time, fall time, width, and AUC. 
    * _TRANSIENTSTREAMLOCS:_ Pre-transient baseline start, pre-transient baseline end, transient peak, rise, and fall locations for each transient to match the cut transient stream data.
    * _TRANSIENTSTREAMDATA: Cut data stream from baseline window start to the end of the post-transient period for each transient event.


**NOTES:**

* If called outside the main _findSessionTransients_ function, note that users will have to input all parameters manually and no defaults will be applied.
* As for the main function, threshold values should be calculated before using the findSessionTransients functions. Typically thresholds are set to 2-3 SDs. If the input data stream is Z scored, this can be the actual SD threshold number. If the input data stream is not Z scored, find the corresponding value to 2-3 SDs for each subject.
* For transient data outputs, each transient is in a separate row. 
* If OUTPUTTRANSIENTDATA is set to anything other than 1, the TRANSIENTSTREAMLOCS and TRANSIENTSTREAMDATA tables will be skipped and not included in the output.


**EXAMPLE:**
```
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients by directly called in the findSessionTransients_blmean function
[data] = findSessionTransients_blmean(data,'sigfiltz_normsession_trimmed','threshold3SD','fs',800,100,2000,0.5,1);
```

### findSessionTransients_blmin
Sub function to detect and quantify transient events in whole session data streams. Transient baselines are defined as the minimum value within the pre-transient baseline window.

This function is called by the main function _findSessionTransients_, which sets default parameters for transient detection and quantification. If called outside the main function, users must specify all input values manually.

**INPUTS:**

* __DATA:__ This is a structure that contains at least the data stream you want to analyze, a field with the threshold values, and a field with the sampling rate of the data stream.
* __WHICHSTREAM:__ The name (string) of the field containing the stream to be analyzed for transients.
* __WHICHTHRESHOLD:__ The name (string) of the field containing the prepared numeric threshold values for each stream. For example, 'threshold_3SD'. 
* __WHICHFS:__ The name (string) of the field containing the sampling rate of the raw data collection in Hz.
* __PREMINSTARTMS:__ Number of millseconds pre-transient to use as the start of the baseline window. _Reccomended value: 800_
* __PREMINENDMS:__ Number of millseconds pre-transient to use as the end of the baseline window. _Reccomended value: 100_
* __POSTTRANSIENTMS:__ Number of millseconds post-transient to use for the post peak baseline and trimmed data output. _Reccomended value: 2000_
* __QUANTIFICATIONHEIGHT:__ The height at which to characterize rise time, fall time, peak width, and AUC. Must be a number between 0 and 1. _Reccomended value: 0.5_
* __OUTPUTTRANSIENTDATA:__ Set to 1 to output cut data streams for each transient event for further analysis or plotting. Set to 0 to skip. _Reccomended value: 1_


**OUTPUTS:**

* __DATA:__ The original data structure with sessiontransients_blmin_THRESHOLDLABEL added in. For more details on individual variables, see the PASTa user guide page _Transient Detection_. The output contains four nested tables: 
    * _INPUTS:_ Includes all function inputs.
    * _TRANSIENTQUANTIFICATION:_ Includes the quantified variables for each transient, including amplitude, rise time, fall time, width, and AUC. 
    * _TRANSIENTSTREAMLOCS:_ Pre-transient baseline, transient peak, rise, and fall locations for each transient to match the cut transient stream data.
    * _TRANSIENTSTREAMDATA: Cut data stream from baseline window start to the end of the post-transient period for each transient event.


**NOTES:**

* If called outside the main _findSessionTransients_ function, note that users will have to input all parameters manually and no defaults will be applied.
* As for the main function, threshold values should be calculated before using the findSessionTransients functions. Typically thresholds are set to 2-3 SDs. If the input data stream is Z scored, this can be the actual SD threshold number. If the input data stream is not Z scored, find the corresponding value to 2-3 SDs for each subject.
* For transient data outputs, each transient is in a separate row. 
* If OUTPUTTRANSIENTDATA is set to anything other than 1, the TRANSIENTSTREAMLOCS and TRANSIENTSTREAMDATA tables will be skipped and not included in the output.


**EXAMPLE:**
```
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients by directly called in the findSessionTransients_blmin function
[data] = findSessionTransients_blmin(data,'sigfiltz_normsession_trimmed','threshold3SD','fs',800,100,2000,0.5,1);
```

### findSessionTransients_localmin
Sub function to detect and quantify transient events in whole session data streams. Transient baselines are defined as the last local minimum value before the transient peak within the pre-transient baseline window.

This function is called by the main function _findSessionTransients_, which sets default parameters for transient detection and quantification. If called outside the main function, users must specify all input values manually.

**INPUTS:**

* __DATA:__ This is a structure that contains at least the data stream you want to analyze, a field with the threshold values, and a field with the sampling rate of the data stream.
* __WHICHSTREAM:__ The name (string) of the field containing the stream to be analyzed for transients.
* __WHICHTHRESHOLD:__ The name (string) of the field containing the prepared numeric threshold values for each stream. For example, 'threshold_3SD'. 
* __WHICHFS:__ The name (string) of the field containing the sampling rate of the raw data collection in Hz.
* __PREMINSTARTMS:__ Number of millseconds pre-transient to use as the start of the baseline window. _Reccomended value: 800_
* __PREMINENDMS:__ Number of millseconds pre-transient to use as the end of the baseline window. _Reccomended value: 100_
* __POSTTRANSIENTMS:__ Number of millseconds post-transient to use for the post peak baseline and trimmed data output. _Reccomended value: 2000_
* __QUANTIFICATIONHEIGHT:__ The height at which to characterize rise time, fall time, peak width, and AUC. Must be a number between 0 and 1. _Reccomended value: 0.5_
* __OUTPUTTRANSIENTDATA:__ Set to 1 to output cut data streams for each transient event for further analysis or plotting. Set to 0 to skip. _Reccomended value: 1_


**OUTPUTS:**

* __DATA:__ The original data structure with sessiontransients_localmin_THRESHOLDLABEL added in. For more details on individual variables, see the PASTa user guide page _Transient Detection_. The output contains four nested tables: 
    * _INPUTS:_ Includes all function inputs.
    * _TRANSIENTQUANTIFICATION:_ Includes the quantified variables for each transient, including amplitude, rise time, fall time, width, and AUC. 
    * _TRANSIENTSTREAMLOCS:_ Pre-transient baseline (local minimum), transient peak, rise, and fall locations for each transient to match the cut transient stream data.
    * _TRANSIENTSTREAMDATA: Cut data stream from baseline window start to the end of the post-transient period for each transient event.


**NOTES:**

* If called outside the main _findSessionTransients_ function, note that users will have to input all parameters manually and no defaults will be applied.
* As for the main function, threshold values should be calculated before using the findSessionTransients functions. Typically thresholds are set to 2-3 SDs. If the input data stream is Z scored, this can be the actual SD threshold number. If the input data stream is not Z scored, find the corresponding value to 2-3 SDs for each subject.
* For transient data outputs, each transient is in a separate row. 
* If OUTPUTTRANSIENTDATA is set to anything other than 1, the TRANSIENTSTREAMLOCS and TRANSIENTSTREAMDATA tables will be skipped and not included in the output.


**EXAMPLE:**
```
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients by directly called in the findSessionTransients_localmin function
[data] = findSessionTransients_localmin(data,'sigfiltz_normsession_trimmed','threshold3SD','fs',800,100,2000,0.5,1);
```

## binSessionTransients
Used to assign individual transient events to time bins within the session for analysis of changes in transients over time. By default, sessions will be divided into 5 minute bins. Number of bins per session will be calculated by the function. Users have the option to override the defaults to change the bin length and/or manually set the number of bins per session.

**INPUTS:**

* __DATA:__ This is a structure that contains at least the data stream you want to analyze, a field with the threshold values, and a field with the sampling rate of the data stream.
* __WHICHSTREAM:__ The name (string) of the field containing the stream to be analyzed for transients.
* __WHICHFS:__ The name (string) of the field containing the sampling rate of the raw data collection in Hz.
* __WHICHTRANSIENTS:__ The name (string) of the parent field containing the table of transietns that you want to identify bins for. This is the output of the _findSessionTransients_ function. The name usually follows the convention *sessiontransients_WHICHBLTYPE_WHICHTHRESHOLD*


**OPTIONAL INPUTS:**
* __WHICHTRANSIENTSTABLE:__ The name (string) of the field within WHICHTRANSIENTS that contains the quantification of individual transient events. This input only needs to be specified if not using the format output from the FINDSESSIONTRANSIENTS functions. _Default: 'transientquantification'_
* __WHICHMAXLOCS:__ NThe name (string) of the field containing the transient max locations (indexes) relative to
the whole session. This input only needs to be specified if not using the format output from the FINDSESSIONTRANSIENTS functions. _Default: 'maxloc'_
* __BINLENGTHMINS:__ Numeric; Bin length in number of minutes. _Default: 5_
* __NBINSOVERRIDE:__ Numeric; Manual override to set the number of bins. If set to anything other than 0, users can override the stream-length based calculation of the number of bins per session and set their own number. _Default: 0_


**OUTPUTS:**

* __DATA:__ The original data structure with sessiontransients_localmin_THRESHOLDLABEL added in. For more details on individual variables, see the PASTa user guide page _Transient Detection_. The output contains four nested tables: 
    * _INPUTS:_ Includes all function inputs.
    * _TRANSIENTQUANTIFICATION:_ Includes the quantified variables for each transient, including amplitude, rise time, fall time, width, and AUC. 
    * _TRANSIENTSTREAMLOCS:_ Pre-transient baseline (local minimum), transient peak, rise, and fall locations for each transient to match the cut transient stream data.
    * _TRANSIENTSTREAMDATA: Cut data stream from baseline window start to the end of the post-transient period for each transient event.


**NOTES:**

* If called outside the main _findSessionTransients_ function, note that users will have to input all parameters manually and no defaults will be applied.
* As for the main function, threshold values should be calculated before using the findSessionTransients functions. Typically thresholds are set to 2-3 SDs. If the input data stream is Z scored, this can be the actual SD threshold number. If the input data stream is not Z scored, find the corresponding value to 2-3 SDs for each subject.
* For transient data outputs, each transient is in a separate row. 
* If OUTPUTTRANSIENTDATA is set to anything other than 1, the TRANSIENTSTREAMLOCS and TRANSIENTSTREAMDATA tables will be skipped and not included in the output.


**EXAMPLE:**
```
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients by directly called in the findSessionTransients_localmin function
[data] = findSessionTransients_localmin(data,'sigfiltz_normsession_trimmed','threshold3SD','fs',800,100,2000,0.5,1);
```

## exportSessionTransients

