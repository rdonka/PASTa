# Function Documentation Overview
This page contains additional documentation for each function within PASTa, as well as examples of inputs.

# Data Preparation Functions
This set of functions is used to prepare raw photometry data, match it with experimental metadata, and load data into a structure in MATLAB. Functions are provided to handle data collected via TDT equipment and software Synapse, or a generic file structure with data streams saved to CSV files.

## loadKeys
Combines subject key and file key into a data structure, and appends the provided computeruserpath to the paths in the file key.

**INPUTS:**

* COMPUTERUSERPATH: A variable containing the unique portion of the file explorer path for the users specific computer. For example, 'C:\Users\rmdon\'. Make sure the computeruserpath ends in a forward slash.
* SUBJECTKEYNAME: A variable containing a string with the name of the subject key csv file for the experiment (see _"3. Data Analysis Pipeline"_ for more details). To omit a subject key and only load in the file key, set subjectkeyname = "".
* FILEKEYNAME: A variable containing a string with the name of the file key csv file for the experiment (see _"3. Data Analysis Pipeline"_ for more details).

**OUTPUTS:**

* EXPERIMENTKEY: A data structure called "experimentkey" that includes the joined file key and subject key with the computer user path appended to raw and extracted folder paths.

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

* RAWFOLDERPATHS: a string array containing the paths to the folder location of the raw data blocks to be extracted. The string array should contain one column with each full path in a separate row.
* EXTRACTEDFOLDERPATHS: a string array containing the paths to the folder location in which to save the extracted MatLab structs for each block to be extracted. The string array should contain one column with each full path in a separate row.  
* SIGSTREAMNAMES: A cell array containing strings with the names of the streams to be treated as signal. Note that only one stream per file can be treated as signal. If different files have different stream names, include all stream names in the cell array.
* BAQSTREAMNAMES: A cell array containing strings with the names of the streams to be treated as background. Note that only one stream per file can be treated as background. If different files have different stream names, include all stream names in the cell array.

**OPTIONAL INPUTS:**

* CLIP: the number of seconds to remove on either end of the data streams. If not specified, defaults to 5 seconds.
* SKIPEXISTING: A binary variable containing a 0 if pre-existing extracted blocks should be re-extracted or a 1 if pre-existing extracted blocks should be skipped. This allows the user to toggle whether or not to extract every block, or only blocks that have not previously been extracted. If not specified, defaults to 1 (skip previously extracted blocks).

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

* EXPERIMENTKEY: A prepared data structure with at minimum the ExtractedFolderPath to locate the individual block structures to be loaded.

**OUTPUTS:**

* DATA: the input data structure with each individual extracted block added by row.

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

* DATA: A data frame containing at least the specified input fields.
* TRIMSTART: The location to start trimming at.
* TRIMEND: The location to end trimming at.
* WHICHSTREAMS: A cell array containing the names of all the streams to be trimmed.

**OPTIONAL INPUTS:**

* WHICHEPOCS: A cell array containing the names of all the epocs to be adjusted due to trimming - subtract the (start loc - 1) from the epoc.

**OUTPUTS:**

* DATA: The data structure with the specified data stream containing the trimmed data.

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

* DATA: A data frame containing at least the specified input fields.
* SIGFIELD: The name (string) of the field containing the signal stream.
* BAQFIELD: The name (string) of the field containing the background stream.
* FS: The sampling rate of the raw data collection in hz.

**OPTIONAL INPUTS:**

* BAQSCALINGTYPE: A string to specify the type of background scaling to apply. Options are 'frequency', 'sigmean', 'OLS', 'detrendOLS', 'smoothedOLS', or 'IRLS'. Default: 'frequency'.
    * 'frequency': Scales the background to the signal channel based on ratio of specified frequency bands in the FFT (frequency domain) of the channels.
    * 'sigmean': Scales the background to the signal channel based on the ratio of the mean of the signal to the mean of the background (time domain).
    * 'OLS': Uses ordinary least-squares regression to generate scaled background.
    * 'detrendOLS': Removes the linear trend from signal and background streams prior to using ordinary least-squares regression to generate scaled background.
    * 'smoothedOLS': Applies lowess smoothing to the signal and background streams prior to using ordinary least-squares regression to generate scaled background.
    * 'IRLS': Uses iteratively reweighted least squares regression to generate scaled background.
* BAQSCALINGFREQ: Only used with 'frequency' scaling. Numeric frequency (Hz) threshold for scaling the background to signal channel. Frequencies above this value will be included in the scaling factor determination. Default: 10 Hz.
* BAQSCALINGPERC: Only used with 'frequency' and 'sigmean' scaling. Adjusts the background scaling factor to be a percednt of the derived scaling factor value. Default: 1 (100%).
* SUBTRACTIONOUTPUT: Output type for the subtracted data. Default: 'dff'
    * 'dff': Outputs subtracted signal as delta F/F.
    * 'df': Outputs subtracted signal as delta F.
* FILTERTYPE: A string to specify the type of filter to apply after subtraction. Default: 'bandpass'.
    * 'nofilter': No filter will be applied.
    * 'bandpass': A bandpass filter will be applied.
    * 'highpass': Only the high pass filter will be applied.
    * 'lowpass': Only the low pass filter will be applied.
* PADDING: Defaults to 1, which applies padding. Padding takes the first 10% of the stream, flips it, and appends it to the data before filtering. Appended data is trimmed after filtration. Set to 0 to turn off padding of data streams. Default: 1.
* PADDINGPERC: Percent of data length to use to determine the number of samples to be appended to the beginning and end of data in padding. Set to minimum 10%. Default: 0.1 (10%).
* FILTERORDER:  The order to be used for the chosen butterworth filter. Default: 3.
* HIGHPASSCUTOFF:  The cutoff frequency (hz) to be used for the high pass butterworth filter. Default: 2.2860.
* LOWPASSCUTOFF: The cutoff to be used for the low pass butterworth filter. Default: 0.0051.
    * NOTE: 'bandpass' applies both the high and low cutoffs to design the filter.
* SUPRESSDISP: If set to anything other than 0, this will suppress the command window displays. Default: 0.

**OUTPUTS:**

* DATA: The original data structure with added fields with the scaled background ('baq_scaled'), subtracted signal ('sigsub'), and subtracted and filtered signal ('sigfilt'). All inputs and defaults will be added to the data structure under the field 'inputs'.
    * NOTE: If using BAQSCALINGMETHOD 'detrendOLS', additional fields containing the detrended signal and background ('sig_detrend' and 'baq_detrend') will be added to the data frame. If using BAQSCALINGMETHOD 'smoothedOLS', additional fields containing the smoothed signal and background ('sig_smoothed' and 'baq_smoothed') will be added to the data frame.


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


## Normalization




# Peak Detection and Quantification Functions
## findpeaks


## peakarea



# Individual Trial Analyses Functions
## cutTrialData

## trialAverage


## trialnormalization

# Plotting Functions
## Whole Session Plots

## FFT Plots

## Peak Plots


## Trial Plots