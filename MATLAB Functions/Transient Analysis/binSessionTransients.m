function [data] = binSessionTransients(data,whichstream,whichfs,whichtransients,varargin)
% BINSESSIONTRANSIENTS      Adds the variable 'Bin' to the transient
%                           quantification table and assigns each transient
%                           a bin based on it's location within the
%                           session.
%
% INPUTS:
%       DATA:               Data structure; This is a structure that contains 
%                           at least the field containing the stream from
%                           which transients were detected, the sampling
%                           rate, and the fields containing the transient
%                           data with a column of transient max indexes.
%
%       WHICHSTREAM:        String; The name of the field containing the 
%                           data stream input to the findSessionTransients
%                           to identify transients from. This is used to
%                           determine how many bins are necessary. For
%                           example, 'sigz_normsession'
%
%       WHICHFS:            String; The name of the field containing the
%                           sampling rate of the streams. For example, 'fs'.
%
%       WHICHTRANSIENTS:    String; The name of the parent field containing 
%                           the table of transients that you want to identify 
%                           bins for. For example, 'sessiontransients_blmin_3SD'.
%
% OPTIONAL INPUTS:
%       WHICHTRANSIENTSTABLE: String; The name of the field within WHICHTRANSIENTS
%                           that contains the quantification of individual 
%                           transient events. This input only needs to be
%                           specified if not using the format output from
%                           the FINDSESSIONTRANSIENTS functions.
%                           Default: 'transientquantification'.
%
%       WHICHMAXLOCS:       String; The name of the field containing the
%                           transient max locations (indexes) relative to
%                           the whole session. This input only needs to be
%                           specified if not using the format output from
%                           the FINDSESSIONTRANSIENTS functions.
%                           Default: 'maxloc'
%
%       BINLENGTHMINS:      Numeric; Bin length in number of minutes.
%                           Default: 5
%
%       NBINSOVERRIDE:      Numeric; Manual override to set the number of 
%                           bins. If set to anything other than 0, 
%                           users can override the stream-length based 
%                           calculation of the number
%                           of bins per session and set their own number.
%                           Default: 0
% OUTPUTS:
%       DATA:               This is the original data structure with bins
%                           added to the specified table of transients.
%
% Written by R M Donka, September 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Settings
% Import required and optional inputs into a structure
    inputs = struct(...
        'whichstream',[],...
        'whichfs',[],...
        'whichtransients',[],...
        'whichtransientstable', [],...
        'whichmaxlocs',[],...
        'binlengthmins',[],...
        'nbinsoverride',[]);
    inputs = parseArgsLite(varargin,inputs);
    
    % Prepare defaults and check for optional inputs
    inputs.whichstream = whichstream;
    inputs.whichfs = whichfs;
    inputs.whichtransients = whichtransients;

    if isempty(inputs.whichtransientstable) % Field containing transient quantification
        whichtransientstable = 'transientquantification';
        inputs.whichtransientstable = whichtransientstable;
    else
        whichtransientstable = inputs.whichtransientstable;
    end

    if isempty(inputs.whichmaxlocs) % Field containing transient max location
        whichmaxlocs = 'maxloc';
        inputs.whichmaxlocs = whichmaxlocs;
    else
        whichmaxlocs = inputs.whichmaxlocs;
    end

    if isempty(inputs.binlengthmins) % Bin length in minutes
        binlengthmins = 5;
        inputs.binlengthmins = binlengthmins;
    else
        binlengthmins = inputs.binlengthmins;
    end

    if isempty(inputs.nbinsoverride) % Bin length in minutes
        nbinsoverride = 0;
        inputs.nbinsoverride = nbinsoverride;
    else
        nbinsoverride = inputs.nbinsoverride;
    end

    disp(append('BIN SESSION TRANSIENTS: Add bin variable to transient quantification table. Binning transients from ', whichtransients, ' into ', num2str(binlengthmins), ' minute bins.')) % Display bin length

    disp('INPUTS:') % Display all input values
    disp(inputs)

    %% Bin Transients
    for eachfile = 1:length(data)
        disp(['   Current File: ',num2str(eachfile)]) % Display which file is being processed

        binsamples = ceil(data(eachfile).(whichfs)*60*binlengthmins); % Pull out n samples per bin

        if nbinsoverride == 0
            nbins = ceil(length(data(eachfile).(whichstream))/binsamples); % Determine number of bins
            disp(append('     ',num2str(nbins),' bins')) % Display which file is being processed
        else
            nbins = nbinsoverride;
            disp(append('     NUMBER OF BINS MANUALLY SET: ', num2str(nbins),' bins')) % Display which file is being processed
        end

        % Add bin transients settings to the transients table
        data(eachfile).(whichtransients).BinSettings.nbins = nbins; % Add number of bins to settings table
        data(eachfile).(whichtransients).BinSettings.nbinsoverride = nbinsoverride; % Add bin number override to settings table
        data(eachfile).(whichtransients).BinSettings.binlengthmins = binlengthmins; % Add bin length in minutes to settings table        
        data(eachfile).(whichtransients).BinSettings.binlengthsamples = binsamples; % Add bin length in samples to settings table

        for eachbin = 0:(nbins-1) % Add bin to transients column
            startbin = (eachbin*binsamples)+1; % Find bin start index
            endbin = (eachbin+1)*binsamples; % Find bin end index
            data(eachfile).(whichtransients).(whichtransientstable).(append('Bin_',num2str(binlengthmins)))((data(eachfile).(whichtransients).(whichtransientstable).(whichmaxlocs) >= startbin & ...
               data(eachfile).(whichtransients).(whichtransientstable).(whichmaxlocs) < endbin)) = eachbin+1; % Add variable for bin to the transient quantification table. Transients are grouped by the bin start and stop indexes. Output column is labeled 'Bin_BINLENGTHMINS'
        end
    end
end