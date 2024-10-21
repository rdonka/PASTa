function [data] = binSessionTransients(data,whichstream,whichfs,whichtransients,whichtransientstable,whichmaxlocs,varargin)
% BINSESSIONTRANSIENTS      Adds the variable 'Bin' to the transient
%                           quantification table
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
%       WHICHTRANSIENTSTABLE: String; The name of the field within WHICHTRANSIENTS
%                           that contains the quantification of individual 
%                           transient events. For example, 'transientquantification'.
%
%       WHICHMAXLOCS:       String; The name of the field containing the
%                           transient max locations (indexes) relative to
%                           the whole session. For example, 'maxloc'.
%
% OPTIONAL INPUTS:
%       BINLENGTHMINS:      Numeric; Bin length in number of minutes.
%                           Default: 5
% OUTPUTS:
%       DATA:               This is the original data structure with bins
%                           added to the specified table of transients
%
% Written by R M Donka, September 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

disp('BIN TRANSIENTS: Add bin variable to transient quantification table.')

%% Prepare Settings
% Import required and optional inputs into a structure
    inputs = struct(...
        'whichstream',[],...
        'whichfs',[],...
        'whichtransients',[],...
        'whichtransientstable', [],...
        'whichmaxlocs',[],...
        'binlengthmins',[]);
    inputs = parseArgsLite(varargin,inputs);
    
    % Prepare defaults and check for optional inputs
    inputs.whichstream = whichstream;
    inputs.whichfs = whichfs;
    inputs.whichtransients = whichtransients;
    inputs.whichtransientstable = whichtransientstable;
    inputs.whichmaxlocs = whichmaxlocs;

    if isempty(inputs.binlengthmins) % Bin length in minutes
        binlengthmins = 5;
        inputs.binlengthmins = binlengthmins;
    else
        binlengthmins = inputs.binlengthmins;
    end

    disp('INPUTS:') % Display all input values
    disp(inputs)

    %% Bin Transients
    for eachfile = 1:length(data)
        disp(['Binning Transients: File ',num2str(eachfile)]) % Display which file is being processed

           data(eachfile).binsamples = floor(data(eachfile).fs*60*minsperbin); % Bin length in samples using sampling rate fs
    
        binsamples = floor(data(eachfile).(whichfs)*60*binlengthmins); % Pull out n samples per bin
        nbins = ceil(length(data(eachfile).(whichstream))/binsamples); % Determine number of bins
        disp(['     Total number of bins: ',num2str(nbins)]) % Display which file is being processed

        for eachbin = 0:(nbins-1) % Add bin to transients column
            startbin = (eachbin*binsamples)+1; % Find bin start index
            endbin = (eachbin+1)*binsamples; % Find bin end index
            data(eachfile).(whichtransients).(whichtransientstable).(append('Bin_',num2str(binlengthmins)))((data(eachfile).(whichtransients).(whichtransientstable).(whichmaxlocs) >= startbin & ...
               data(eachfile).(whichtransients).(whichtransientstable).(whichmaxlocs) < endbin)) = eachbin+1; % Add variable for bin to the transient quantification table. Transients are grouped by the bin start and stop indexes. Output column is labeled 'Bin_BINLENGTHMINS'
        end
    end
end