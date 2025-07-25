function defaultparameters = configDefaultParameters(callerFunction)
% CONFIGDEFAULTPARAMETERS  Returns a struct of default parameters for various functions.
%
%   DEFAULTPARAMETERS = CONFIGDEFAULTPARAMETERS(CALLERFUNCTION) returns only
%   the parameters relevant to the specified function. If no function name is
%   provided, all parameters are returned.
%
% INPUTS:
%   CALLERFUNCTION - (optional) String; Name of the function requesting parameters.
%                    Example: 'cropFPdata', 'extractTDTdata'
%
% OUTPUTS:
%   DEFAULTPARAMETERS - Struct containing relevant parameters.
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

% Define all default parameters for different functions
    % extractTDTdata
    allparameters.extractTDTdata.trim = 5;
    allparameters.extractTDTdata.skipexisting = 1;

    % extractCSVdata
    allparameters.extractCSVdata.trim = 5;
    allparameters.extractCSVdata.skipexisting = 1;
    allparameters.extractCSVdata.loadepocs = 0;

    % subtractFPdata
    allparameters.subtractFPdata.baqscalingtype = 'frequency'; % Subtraction type defaults to frequency scaling
    allparameters.subtractFPdata.baqscalingfreqmin = 10; % Background scaling frequency min defaults to 10 Hz
    allparameters.subtractFPdata.baqscalingfreqmax = 100; % Background scaling frequency max defaults to 100 Hz
    allparameters.subtractFPdata.baqscalingperc = 1; % Background scaling percent defaults to 100%
    allparameters.subtractFPdata.subtractionoutput = 'dff';  % Subtraction output defaults to delta F/F
    allparameters.subtractFPdata.artifactremoval = 0; % Artifact removal defaults to false
    allparameters.subtractFPdata.filtertype = 'bandpass'; % Filter type defaults to bandpass filter
    allparameters.subtractFPdata.padding = 1; % Pre-filter padding length defaults to 10% of stream length
    allparameters.subtractFPdata.paddingperc = 0.1; % Filter type defaults to bandpass filter
    allparameters.subtractFPdata.filterorder = 3; % Filter order defaults to 3rd order
    allparameters.subtractFPdata.highpasscutoff = 0.0051; % High pass frequency cuttoff defaults to 0.0051 Hz
    allparameters.subtractFPdata.lowpasscutoff = 2.2860; % Low pass frequency cuttoff defaults to 2.2860 Hz

    % removeStreamArtifacts
    allparameters.removeStreamArtifacts.outlierthresholdk = 10; % Outlier detection threshold; k = 3 or 4 is a common choice for outlier detection
    allparameters.removeStreamArtifacts.artifactremovalwindow = 0.3; % Seconds to replace with NaNs before and after artifact
    allparameters.removeStreamArtifacts.artifactampthreshold_max = 8; % Artifact SD detection threshold - max artifacts
    allparameters.removeStreamArtifacts.artifactampthreshold_min = 8; % Artifact SD detection threshold - min artifacts
    allparameters.removeStreamArtifacts.bucketsizeSecs = 30; % Bucket size (seconds) for session amplitude and mean calculations

    % findSessionTransients
    allparameters.findTransients.bltype = 'blmean'; % Default to baseline mean method
    allparameters.findTransients.blstartms = 1000; % Pre-transient baseline ms start
    allparameters.findTransients.blendms = 200; % Pre-transient baseline ms end
    allparameters.findTransients.posttransientms = 2000; % Post-transient fall window
    allparameters.findTransients.AUCwindowms = 1000; % Window size (ms) for AUC window
    allparameters.findTransients.compoundtransientwindowms = 2000; % Window size to search before and after each event for compound transients
    allparameters.findTransients.quantificationheight = 0.5; % Height for quantification of transients (rise/fall/AUC)
    allparameters.findTransients.outputtransientdata = 1; % Logical; If set to 1 (true), data streams for individual transients will be added to the data structure
    allparameters.findTransients.outputpremaxS = 3; % Logical; If set to 1 (true), data streams for individual transients will be added to the data structure
    allparameters.findTransients.outputpostmaxS = 5; % Logical; If set to 1 (true), data streams for individual transients will be added to the data structure

    % binSessionTransients
    allparameters.binTransients.binlengthmins = 5; % Bin length (mins) for transient analysis
    allparameters.binTransients.nbinsoverride = 0; % Mannual override to set number of bins; Set to 0 to calculate number of bins based on length of session    
    allparameters.binTransients.manuallydefinebins = false; % Mannually specify bin start and end indexes

    % exportSessionTransients
    allparameters.exportTransients.transientquantificationfieldname = 'transientquantification'; % Name of field containing table of quantified transients
    allparameters.exportTransients.exportfilename = ''; % Empty file name to trigger automatic naming

    % plotTraces
    allparameters.plotTraces.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotTraces.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotTraces.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotNormTraces
    allparameters.plotNormTraces.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotNormTraces.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotNormTraces.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotFFTpower
    allparameters.plotFFTpower.xmax = 100; % X axis cutoff (hz) for plots
    allparameters.plotFFTpower.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotFFTpower.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotFFTpower.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotFFTmag
    allparameters.plotFFTmag.xmax = 100; % X axis cutoff (hz) for plots
    allparameters.plotFFTmag.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotFFTmag.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotFFTmag.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotTransients
    allparameters.plotTransients.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotTransients.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotTransients.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotTransientBins
    allparameters.plotTransientBins.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotTransientBins.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotTransientBins.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotTransientTraces
    allparameters.plotTransientTraces.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotTransientTraces.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotTransientTraces.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

    % plotTransientTraceBins
    allparameters.plotTransientTraceBins.saveoutput = 0; % Logical; Set to 1 to automatically save plot to plotfilepath
    allparameters.plotTransientTraceBins.outputfiletype = 'png'; % String; If saveoutput = 1, save plot as png
    allparameters.plotTransientTraceBins.plotfilepath = ''; % Empty file name - replace with manual input if saveoutput set to 1

% If a specific function is requested, return only its parameters
    if nargin > 0 && isfield(allparameters, callerFunction)
        defaultparameters = allparameters.(callerFunction);
    else
        % If no function is specified, return the entire structure
        defaultparameters = allparameters;
    end
end
