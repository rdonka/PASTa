function [data] = findSessionTransients(data,whichbltype,whichstream,whichthreshold,whichfs,varargin)
% FINDSESSIONTRANSIENTS_BLMIN   Finds transients for the whole session. Use
%                               of this main function will set default values
%                               for parameters not specific as optional inputs.
%
%                               For additional details for each baseline
%                               method, see subfunctions:
%                               FINDSESSIONTRANSIENTS_BLMIN
%                               FINDSESSIONTRANSIENTS_BLMEAN
%                               FINDSESSIONTRANSIENTS_LOCALMIN
%                               
% INPUTS:
%       DATA:           This is a structure that contains at least the data
%                       stream you want to analyze. For example, 
%                       'sigfiltz_normsession'.
%
%
%       WHICHBLTYPE:    A variable containing a string with the type of
%                       pre-transient baseline to use for amplitude
%                       inclusion and quantification.
%                       OPTIONS:
%                           'blmin': Pre-transient baselines are set to the
%                               minimum value within the pre-transient window.
%                           'blmean': Pre-transient baselines are set to the
%                               mean of the pre-transient window.
%                           'localmin': Pre-transient baselines are set to
%                               the local minimum directly preceding the
%                               transient within the baseline window.
%                       Default: 'blmin'
%
%       WHICHSTREAM:    A variable containing a string with the name of the 
%                       field containing the stream to be analyzed for 
%                       transients. For example, 'sigfiltz_normsession'.
%
%       WHICHTHRESHOLD: A variable containing a string with the name of the
%                       field containing the prepared numeric threshold
%                       values for each stream. For example, 'threshold_3SD'.
%                       NOTE: Threshold values should be calculated
%                       before using the findSessionTransients functions.
%                       Typically thresholds are set to 2-3 SDs. If the
%                       input data stream is Z scored, this can be the
%                       actual SD threshold number. If the input data
%                       stream is not Z scored, find the corresponding
%                       value to 2-3 SDs for each subject.
%                     
%       WHICHFS:        The name of the field containing the sampling rate
%                       of the streams (fs).
%
%
% OPTIONAL INPUTS:
%       PREMINSTARTMS:  Number of millseconds pre-transient to use as the
%                       start of the baseline window.
%                       Default: 1000
%
%       PREMINENDMS:    Number of millseconds pre-transient to use as the
%                       end of the baseline window.
%                       Default: 100
%
%       POSTTRANSIENTMS: Number of millseconds post-transient to use for
%                       the post peak baseline and trimmed data output.
%                       Default: 2000
%
%       QUANTIFICATIONHEIGHT: The height at which to characterize rise time,
%                       fall time, peak width, and AUC. Must be a number 
%                       between 0 and 1. Default: 0.5
%  
%       OUTPUTTRANSIENTDATA: Set to 1 to output cut data streams for each
%                       transient event. Set to 0 to skip.
%                       Default: 1
%
% OUTPUTS:
%       DATA:           The original data structure with
%                       sessiontransients_WHICHBLTYPE_THRESHOLDLABEL added in.
%                       For more details on individual variables, see the
%                       PASTa user guide. 
%                       The output contains four nested tables: 
%                       INPUTS: Includes all required and optional inputs.
%                           If optional inputs are not specified, defaults
%                           will be applied.
%                       TRANSIENTQUANTIFICATION: Includes the quantified
%                           variables for each transient, including
%                           amplitude, rise time, fall time, width, and
%                           AUC. 
%                       TRANSIENTSTREAMLOCS: Pre-transient baseline, 
%                           transient peak, rise, and fall locations for 
%                           each transient to match the cut transient 
%                           stream data.
%                       TRANSIENTSTREAMDATA: Cut data stream from baseline
%                           start to the end of the post-transient period
%                           for each transient event.
%                       Note that for all data outputs, each transient is
%                       in a separate row. If OUTPUTTRANSIENTDATA is set to 
%                       anything other than 1, the TRANSIENTSTREAMLOCS and
%                       TRANSIENTSTREAMDATA tables will be skipped and not
%                       included in the output.
%
% Written by R M Donka, October 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Settings
% Import required and optional inputs into a structure
    maininputs = struct(...
        'whichbltype',[],...
        'whichstream',[],...
        'whichthreshold',[],...
        'whichfs',[],....
        'preminstartms', [],...
        'preminendms',[],...
        'posttransientms',[],...
        'quantificationheight',[],...
        'outputtransientdata',[]);
    maininputs = parseArgsLite(varargin,maininputs);
    
    % Prepare defaults and check for optional inputs
    maininputs.whichbltype = whichbltype;
    maininputs.whichstream = whichstream;
    maininputs.whichthreshold = whichthreshold;
    maininputs.whichfs = whichfs;

    if isempty(maininputs.preminstartms) % Pre-transient baseline ms start
        preminstartms = 1000;
        maininputs.preminstartms = preminstartms;
    else
        preminstartms = maininputs.preminstartms;
    end
    if isempty(maininputs.preminendms) % Pre-transient baseline ms end
        preminendms = 100;
        maininputs.preminendms = preminendms;
    else
        preminendms = maininputs.preminendms;
    end
    if isempty(maininputs.posttransientms) % Post transient ms for fall/AUC and output data
        posttransientms = 2000;
        maininputs.posttransientms = posttransientms;
    else
        posttransientms = maininputs.posttransientms;
    end
    if isempty(maininputs.quantificationheight) % Height for rise/fall/AUC; defaults to 0.5 (half height)
        quantificationheight = 0.5; 
        maininputs.quantificationheight = quantificationheight;
    else
        quantificationheight = maininputs.quantificationheight;
    end
    if isempty(maininputs.outputtransientdata)
        outputtransientdata = 1; % If set to 1, data streams for individual transients will be added to the data structure
        maininputs.outputtransientdata = outputtransientdata;
    else
        outputtransientdata = maininputs.outputtransientdata;
    end

%% Call findSessionTransients subfunctions: findSessionTransients_blmin, findSessionTransients_blmean, or findSessionTransients_localmin
    if strcmp(whichbltype,'blmin') == true
        [data] = findSessionTransients_blmin(data,whichstream,whichthreshold,whichfs,preminstartms,preminendms,posttransientms,quantificationheight,outputtransientdata);
    elseif strcmp(whichbltype,'blmean') == true
        [data] = findSessionTransients_blmean(data,whichstream,whichthreshold,whichfs,preminstartms,preminendms,posttransientms,quantificationheight,outputtransientdata);
    elseif strcmp(whichbltype,'localmin') == true
        [data] = findSessionTransients_localmin(data,whichstream,whichthreshold,whichfs,preminstartms,preminendms,posttransientms,quantificationheight,outputtransientdata);
    else
        disp("ERROR: No viable transient baseline type specific. WHICHBLTYPE must be set to 'blmin', 'blmean', or 'localmin'")
    end
end