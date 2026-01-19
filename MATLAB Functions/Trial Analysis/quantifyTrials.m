function [data] = quantifyTrials(data,trialstreamfieldname,trialphasestartfieldname,trialphaseendfieldname,fsfieldname,varargin)
% QUANTIFYTRIALS  Quantifies trial data streams by trial events provided in
%                 'trialphasestartfieldname' and 'trialphaseendfieldname'.
%                 Event epochs must be provided relative to the trial data
%                 streams. Multiple windows within a trial can be provided
%                 by passing arrays of epochs to 'trialphasestartfieldname' 
%                 and 'trialphaseendfieldname'.
% 
%
%   QUANTIFYTRIALS(DATA, TRIALSTREAMFIELDNAME, TRIALPHASESTARTFIELDNAME, TRIALPHASEENDFIELDNAME,...
%                   FSFIELDNAME,'PARAM1', VAL1, ...)
%   quantifies the mean, max, and AUC of individual trial streams and adds
%   the field 'trialquantification' to the data frame.
%
% REQUIRED INPUTS:
%   DATA                        - Struct array; Must contain at least a field with the desired cut
%                                 trial streams to be quantified with each trial as a row.
%
%   TRIALSTREAMFIELDNAME        - String; The name of the field in the data structure that
%                                 contains the trial streams to be quantified.
%
%   TRIALPHASESTARTFIELDNAME    - String; The name of the field containing epochs to provide the 
%                                 start of each phase window for quantification (i.e., baseline, 
%                                 stimulus delivery, post stimulus delivery).
%
%   TRIALPHASEENDFIELDNAME     - String; The name of the field containing epochs to provide the 
%                                end of each phase window for quantification.
%
%   FSFIELDNAME                 - String; name of the field in DATA containing the
%                                 sampling rate (e.g., 'fs').
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%   'trialphaselabels'    - Cell array of strings; names to label each phase provided by 
%                           the trial event start and end epochs.
%
%   'trialidfieldnames'   - String; name of field containing trial ID variable (i.e., trial type) 
%                           to be added to trial quantification table.
%
% OUTPUTS:
%   DATA                - Original data structure with a field added
%                         containing the trial quantification. The output
%                         field name will be 'trialquantification_<TRIALSTREAMFIELDNAME>'
%
% EXAMPLE:
%   data = quantifyTrials(data, 'trial_sigfiltz_normbaseline', 'trial_phasestartidxs', 'trial_phasendidxs', 'fs');
%
% See also: cutTrialdata, centerTrialdata, normTrialdata,
% exportTrialquantification
%
% Author:  Rachel Donka (2026)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/


    %% Prepare Settings
    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'trialphaselabels', []); % exportfilepath: string with the custom exportfilepath
    addParameter(p, 'trialidfieldnames', []); 

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Initialize params
    params.trialstreamfieldname = trialstreamfieldname;
    params.trialphasestartfieldname = trialphasestartfieldname;
    params.trialphaseendfieldname = trialphaseendfieldname;
    
    if ~isempty(params.trialphaselabels)
        trialphaselabels = params.trialphaselabels;
    end

    if ~isempty(params.trialidfieldnames)
        trialidfieldnames = params.trialidfieldnames;
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
            
            % Prepare phase epochs
            trialeventstartepochs = data(eachfile).(trialphasestartfieldname)(eachtrial,:);
            trialeventendepochs = data(eachfile).(trialphaseendfieldname)(eachtrial,:);

            if length(trialeventstartepochs) ~= length(trialeventendepochs)
                error('ERROR: Number of input trial start epochs not the same as number of input trial end epochs.')
            end

            % Quantify trial stream for each phase
            for eachphase = 1:length(trialeventstartepochs)
                currvals = [];
                currtrialeventstartidx = trialeventstartepochs(eachphase);
                currtrialeventendidx = trialeventendepochs(eachphase);

                currtrialphasedata = data(eachfile).(trialstreamfieldname)(eachtrial,currtrialeventstartidx:currtrialeventendidx);

                % Find duration of phase
                currtrialphasesamples = size(currtrialphasedata,2);
                currtrialphaseS = currtrialphasesamples/fs;

                % Prepare variables
                currvals.Trial = eachtrial;
                currvals.Stream = {trialstreamfieldname};
                currvals.PhaseNum = eachphase;

                % Add phase label if optional input 'trialphaselabels' provided
                if ~isempty(params.trialphaselabels)
                    currvals.Phase = {trialphaselabels{eachphase}};
                end

                % Add trial phase samples
                currvals.TrialPhaseS = currtrialphaseS;

                % Add trial phase quantification
                currvals.TrialPhasemean = mean(currtrialphasedata,2,'omitnan'); % Phase mean
                currvals.TrialPhasemin = min(currtrialphasedata,[],2,'omitnan'); % Phase min
                currvals.TrialPhasemax = max(currtrialphasedata,[],2,'omitnan'); % Phase max
                currvals.TrialPhaserange = currvals.TrialPhasemax - currvals.TrialPhasemin; % Phase range
                currvals.TrialPhasesd = std(currtrialphasedata,'omitnan'); % Phase standard deviation

                % Add amplitude relative to baseline
                if eachphase == 1 % Skip for baseline phase
                    currvals.TrialPhaseMaxAmp = nan;
                    currvals.TrialPhaseMeanAmp = nan; 
                    currvals.TrialPhaseAUC = nan;
                else % Determine amplitude for all phases after baseline
                    currtrialphaseAUCdata = currtrialphasedata - currtrialvals.TrialPhasemean(1); 

                    currvals.TrialPhaseMaxAmp = currvals.TrialPhasemax - currtrialvals.TrialPhasemax(1); % Phase max amplitude relative to baseline
                    currvals.TrialPhaseMeanAmp = currvals.TrialPhasemean - currtrialvals.TrialPhasemean(1); % Phase mean amplitude relative to baseline
                    currvals.TrialPhaseAUC = trapz(currtrialphaseAUCdata, 2)/currtrialphasesamples; % Phase AUC relative to baseline
                end

                % Add trial ids if optional input 'trialidfieldnames' provided
                if ~isempty(params.trialidfieldnames)
                    for eachtrialidfieldnames = 1:length(trialidfieldnames)
                        currtrialidfieldname = trialidfieldnames{eachtrialidfieldnames};
                        currvals.(currtrialidfieldname) = data(eachfile).(currtrialidfieldname)(eachtrial);
                    end
                end
                currvalstable = struct2table(currvals);
                currtrialvals = [currtrialvals; currvalstable];
            end
            currfiletrialvals = [currfiletrialvals; currtrialvals];
        end
    data(eachfile).(['trialquantification_',trialstreamfieldname]) = currfiletrialvals; % Add trial quantification to data structure
    end
end


% Copyright (C) 2026 Rachel Donka
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