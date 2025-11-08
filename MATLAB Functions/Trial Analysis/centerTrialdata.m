function [data] = centerTrialdata(data,trialdatafieldname,trialBLstartepocfieldname,trialBLendepocfieldname)
% CUTTRIALDATA    Cuts fiber photometry signal data into individual trials based on epocs.
%                 Requires a data frame containing at least the specified signal and relevant
%                 epocs for each trial.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       desired signal to be analyzed and epocs to mark each
%                       trial.
%
%       trialdatafieldname:    The name of the field in the data structure that
%                       contains the signal to be cut.
%
%       trialBLstartepocfieldname:      The name of the epoc to provide the start index for
%                       each trial.
%
%       trialBLendepocfieldname:        The name of the epoc to provide the end index for
%                       each trial.
%
%       epocsfieldnames:      A list of the names of the fields of any other epocs
%                       relevant to trial analysis.
% OUTPUTS:
%       DATA:           This is the original data structure with added
%                       fields for trial signal data and epocs.
%
%                       data.trial_streamfieldname: contains the cut trial data
%                       with each trial saved as a row.
%
%                       data.trial_startepocfieldname: contains conthe index for the epoc
%                       specified as epoc start relative to the trimmed
%                       trial data.
%
%                       data.trial_whichepocs: each other epoc will generate
%                       a field that contains the index for that epoc
%                       relative to the startepocfieldname for each trial.
%
% Written by R M Donka, January 2024
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

% Set up optional inputs
%optionalinputs = struct('epocsfieldnames',[]);
%optionalinputs = parseArgsLite(varargin,optionalinputs);
    

% Cut trial data
    disp(['CENTERING TRIAL DATA: Centering ', trialdatafieldname, ' to trial baseline.']);
    disp(['   Baseline start indices field: ', trialBLstartepocfieldname]);
    disp(['   Baseline end indices field: ', trialBLendepocfieldname]);

    for eachfile = 1:length(data)
        ntrials = size(data(eachfile).(trialdatafieldname),1);
        maxtrialsamples = length(data(eachfile).(trialdatafieldname));
        % Prepare data structure
        data(eachfile).(append(trialdatafieldname,'centered'))(1:ntrials,1:maxtrialsamples) = NaN;

        for eachtrial = 1:ntrials
            currtrialsamples = length(data(eachfile).(trialdatafieldname)(eachtrial,:));
            currtrialBLmean = mean(data(eachfile).(trialdatafieldname)(eachtrial,data(eachfile).(trialBLstartepocfieldname)(eachtrial):data(eachfile).(trialBLendepocfieldname)(eachtrial)),'omitnan');

            data(eachfile).(append(trialdatafieldname,'centered'))(eachtrial,1:currtrialsamples) = data(eachfile).(trialdatafieldname)(eachtrial,:) - currtrialBLmean; % Z score the trial to the trial baseline
        end
    end
end