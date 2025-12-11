function [data] = quantifyTrials(data,trialstreamfieldname,trialeventstartfieldname,trialeventendfieldname,fsfieldname,varargin)
% QUANTIFYTRIALS  Assigns each transient to a time bin within the session.
%
%   QUANTIFYTRIALS(DATA, TRIALSTREAMFIELDNAME, 'PARAM1', VAL1, ...)
%   quantifies the mean, max, and AUC of individual trial streams and adds
%   the field 'trialquantification' to the data frame.
%
% REQUIRED INPUTS:
%       DATA       - Data structure containing at least a field with cut
%                    streams of trial data with each trial as a row.
%
%       TRIALSTREAMFIELDNAME    - String; name of the field in DATA containing
%                                 the cut trial data to be quantified. Each
%                                 trial should be in each row.
%
%       TRIALEVENTSTARTFIELDNAME - String; name of the field in DATA
%                                  containing the phase start epochs spatially 
%                                  matched to the cut trial data specified by 
%                                  TRIALSTREAMFIELDNAME.
%
%       TRIALEVENTENDFIELDNAME - String; name of the field in DATA
%                                  containing the phase end epochs spatially 
%                                  matched to the cut trial data specified by 
%                                  TRIALSTREAMFIELDNAME.
%
%       FSFIELDNAME         - String; name of the field in DATA containing the sampling rate (fs)
%                             of the data stream.
%
% OPTIONAL INPUTS
%       TRIALPHASELABELS    - Cell array of strings; labels for each phase
%                             of the trial structure.
% OUTPUTS:
%   DATA                - Original data structure with a field added
%                         containing the trial quantification. The output
%                         field name will be 'trialquantification_<TRIALSTREAMFIELDNAME>'
%
% EXAMPLE:
%   data = quantifyTrials(data, 'trial_sigfiltz_normbaseline', 'trial_solO', 'trial_solF', 'fs');
%
% See also: findTransients
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/


    %% Prepare Settings
    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'trialphaselabels', []); % exportfilepath: string with the custom exportfilepath

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Initialize params
    params.trialstreamfieldname = trialstreamfieldname;
    params.trialeventstartfieldname = trialeventstartfieldname;
    params.trialeventendfieldname = trialeventendfieldname;
    
    if ~isempty(params.trialphaselabels)
        trialphaselabels = params.trialphaselabels;
    end

    % Display
    disp(['QUANTIFYTRIALS: Add trial quantification to data structure. Trial streams will be quantified from the field: ', trialstreamfieldname]) % Display bin length
    disp('   PARAMETERS:') % Display all parameters
    disp(params)


    %% QUANTIFY TRIAL RESPONSES
    for eachfile = 1:length(data)
        disp(['     QUANTIFYING TRIALS: File ', num2str(eachfile)])
        currfiletrialvals = table();
    
        ntrials = size(data(eachfile).(trialstreamfieldname),1);        
        fs = data(eachfile).(fsfieldname);
    
        for eachtrial = 1:ntrials
            % Prepare table
            currtrialvals = table();
            
            trialeventstartepochs = data(eachfile).(trialeventstartfieldname)(eachtrial,:);
            trialeventendepochs = data(eachfile).(trialeventendfieldname)(eachtrial,:);

            if length(trialeventstartepochs) ~= length(trialeventendepochs)
                error('ERROR: Number of input trial start epochs not the same as number of input trial end epochs.')
            end

            for eachphase = 1:length(trialeventstartepochs)
                currvals = [];
                currtrialeventstartidx = trialeventstartepochs(eachphase);
                currtrialeventendidx = trialeventendepochs(eachphase);

                currtrialphasedata = data(eachfile).(trialstreamfieldname)(eachtrial,currtrialeventstartidx:currtrialeventendidx);

                % Find duration of phase
                currtrialphasesamples = size(currtrialphasedata,2);
                currtrialphaseS = currtrialphasesamples/fs;

                currvals.Trial = eachtrial;
                currvals.Stream = {trialstreamfieldname};
                currvals.PhaseNum = eachphase;
                if ~isempty(params.trialphaselabels)
                    currvals.Phase = {trialphaselabels{eachphase}};
                end
                currvals.TrialPhaseS = currtrialphaseS;

                currvals.TrialPhasemean = mean(currtrialphasedata,2,'omitnan');
                currvals.TrialPhasemin = min(currtrialphasedata,[],2,'omitnan');
                currvals.TrialPhasemax = max(currtrialphasedata,[],2,'omitnan');
                currvals.TrialPhaserange = currvals.TrialPhasemax - currvals.TrialPhasemin;
                currvals.TrialPhasesd = std(currtrialphasedata,'omitnan');

                if eachphase == 1
                    currvals.TrialPhaseMaxAmp = nan;
                    currvals.TrialPhaseMeanAmp = nan;
                    currvals.TrialPhaseAUC = nan;
                else
                    currtrialphaseAUCdata = currtrialphasedata - currtrialvals.TrialPhasemean(1);

                    currvals.TrialPhaseMaxAmp = currvals.TrialPhasemax - currtrialvals.TrialPhasemax(1);
                    currvals.TrialPhaseMeanAmp = currvals.TrialPhasemean - currtrialvals.TrialPhasemean(1);
                    currvals.TrialPhaseAUC = trapz(currtrialphaseAUCdata, 2)/currtrialphasesamples;
                end
                currvalstable = struct2table(currvals);
                currtrialvals = [currtrialvals; currvalstable];
            end
            currfiletrialvals = [currfiletrialvals; currtrialvals];
        end
    data(eachfile).(['trialquantification_',trialstreamfieldname]) = currfiletrialvals;
    end
end


% Copyright (C) 2025 Rachel Donka
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.