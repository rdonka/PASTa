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
%                                  containing the event start epochs spatially 
%                                  matched to the cut trial data specified by 
%                                  TRIALSTREAMFIELDNAME.
%
%       TRIALEVENTENDFIELDNAME - String; name of the field in DATA
%                                  containing the event end epochs spatially 
%                                  matched to the cut trial data specified by 
%                                  TRIALSTREAMFIELDNAME.
%
%       FSFIELDNAME         - String; name of the field in DATA containing the sampling rate (fs)
%                             of the data stream.
%
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
    % Initialize params
    params = struct();
    params.trialstreamfieldname = trialstreamfieldname;
    params.trialeventstartfieldname = trialeventstartfieldname;
    params.trialeventendfieldname = trialeventendfieldname;
    
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
            currvals = table();
    
            % Subset data for pre event, event, and post event periods
            currtrialpreeventdata = data(eachfile).(trialstreamfieldname)(eachtrial,1:data(eachfile).(trialeventstartfieldname)(eachtrial)-1);
            currtrialeventdata = data(eachfile).(trialstreamfieldname)(eachtrial,data(eachfile).(trialeventstartfieldname)(eachtrial):data(eachfile).(trialeventendfieldname)(eachtrial));
            currtrialposteventdata = data(eachfile).(trialstreamfieldname)(eachtrial,data(eachfile).(trialeventendfieldname)(eachtrial):end);
            currtrialeventAUCdata = currtrialeventdata - mean(currtrialpreeventdata,'omitnan');
                
            % Find duration of each period
            currtrialpreeventsamples = size(currtrialpreeventdata,2);
            currtrialpreeventS = currtrialpreeventsamples/fs;
    
            currtrialeventsamples = size(currtrialeventdata,2);
            currtrialeventS = currtrialeventsamples/fs;
    
            currtrialposteventsamples = size(currtrialposteventdata,2);
            currtrialposteventS = currtrialposteventsamples/fs;
    
            currvals.Trial = eachtrial;
            currvals.Stream = {trialstreamfieldname};
            currvals.TrialPreEventS = currtrialpreeventS;
            currvals.TrialEventS = currtrialeventS;
            currvals.TrialPostEventS = currtrialposteventS;
    
            currvals.TrialPreEventmean = mean(currtrialpreeventdata,2,'omitnan');
            currvals.TrialPreEventmin = min(currtrialpreeventdata,[],2,'omitnan');
            currvals.TrialPreEventmax = max(currtrialpreeventdata,[],2,'omitnan');
            currvals.TrialPreEventrange = currvals.TrialPreEventmax - currvals.TrialPreEventmin;
            currvals.TrialPreEventsd = std(currtrialpreeventdata,'omitnan');
    
            currvals.TrialEventmean = mean(currtrialeventdata,2,'omitnan');
            currvals.TrialEventmin = min(currtrialeventdata,[],2,'omitnan');
            currvals.TrialEventmax = max(currtrialeventdata,[],2,'omitnan');
            currvals.TrialEventrange = currvals.TrialEventmax - currvals.TrialEventmin;
            currvals.TrialEventsd = std(currtrialeventdata,'omitnan');
            currvals.TrialEventMaxAmp = currvals.TrialEventmax - currvals.TrialPreEventmean;
            currvals.TrialEventMeanAmp = currvals.TrialEventmean - currvals.TrialPreEventmean;
            currvals.TrialEventAUC = trapz(currtrialeventAUCdata, 2);
    
            currvals.TrialPostEventmean = mean(currtrialposteventdata,2,'omitnan');
            currvals.TrialPostEventmin = min(currtrialposteventdata,[],2,'omitnan');
            currvals.TrialPostEventmax = max(currtrialposteventdata,[],2,'omitnan');
            currvals.TrialPostEventrange = currvals.TrialPostEventmax - currvals.TrialPostEventmin;
            currvals.TrialPostEventsd = std(currtrialposteventdata,'omitnan');
            currvals.TrialPostEventMaxAmp = currvals.TrialPostEventmax - currvals.TrialPreEventmean;
            currvals.TrialPostEventMeanAmp = currvals.TrialPostEventmean - currvals.TrialPreEventmean;
    
            currvals.Trial = eachtrial;
            currvals.Stream = {trialstreamfieldname};
            currvals.TrialPreEventS = currtrialpreeventS;
            currvals.TrialEventS = currtrialeventS;
            currvals.TrialPostEventS = currtrialposteventS;
    
            currfiletrialvals = [currfiletrialvals; currvals];
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