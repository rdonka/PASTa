function [data] = normTrialdata(data,trialstreamfieldname,trialBLstartepocfieldname,trialBLendepocfieldname)
% NORMTRIALDATA   Normalizing individual trial data streams prepared with cutTrialdata to trial baseline period.
%                 Requires a data frame containing at least the specified
%                 trial data streams and trial baseline start and end event
%                 epochs for each trial.
%
% REQUIRED INPUTS:
%   DATA                        - Struct array; Must contain at least the desired cut trial stream 
%                                 to be normalized and epocs to mark the start and end of each trial's
%                                 baseline period.
%
%   TRIALSTREAMFIELDNAME          - String; The name of the field in the data structure that
%                                 contains the signal to be cut.
%
%   TRIALBLSTARTEPOCFIELDNAME   - String; The name of the epoc to provide the start index for
%                                 the baseline period of each trial.
%
%   TRIALBLENDEPOCFIELDNAME     - String; The name of the epoc to provide the end index for
%                                 the baseline period of each trial.
%
% OUTPUT:
%   DATA    - Struct array; This is the original data structure with an added
%             field for the normalized trial signal data.
%               - data.trial_<trialstreamfieldname>z: contains the normalized
%                 (z score) trial data stream.
%
% Author:  Rachel Donka (2026)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/


% Cut trial data
    disp(['NORMALIZING TRIAL DATA: Normalizing ', trialstreamfieldname, ' to trial baseline.']);
    disp(['   Baseline start indices field: ', trialBLstartepocfieldname]);
    disp(['   Baseline end indices field: ', trialBLendepocfieldname]);

    for eachfile = 1:length(data)
        ntrials = size(data(eachfile).(trialstreamfieldname),1);
        maxtrialsamples = length(data(eachfile).(trialstreamfieldname));
        % Prepare data structure
        data(eachfile).(append(trialstreamfieldname,'z'))(1:ntrials,1:maxtrialsamples) = NaN;

        % Normalize trials
        for eachtrial = 1:ntrials
            currtrialsamples = length(data(eachfile).(trialstreamfieldname)(eachtrial,:));
            currtrialBLmean = mean(data(eachfile).(trialstreamfieldname)(eachtrial,data(eachfile).(trialBLstartepocfieldname)(eachtrial):data(eachfile).(trialBLendepocfieldname)(eachtrial)),'omitnan');
            currtrialBLsd = std(data(eachfile).(trialstreamfieldname)(eachtrial,data(eachfile).(trialBLstartepocfieldname)(eachtrial):data(eachfile).(trialBLendepocfieldname)(eachtrial)),'omitnan');

            data(eachfile).(append(trialstreamfieldname,'z'))(eachtrial,1:currtrialsamples) = (data(eachfile).(trialstreamfieldname)(eachtrial,:) - currtrialBLmean)/currtrialBLsd; % Z score the trial to the trial baseline
        end
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