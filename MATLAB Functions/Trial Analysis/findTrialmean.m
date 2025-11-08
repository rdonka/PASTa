function [data] = findTrialMean(data,trialstreamfieldname,varargin)
% CUTTRIALDATA    Cuts fiber photometry signal data into individual trials based on epocs.
%                 Requires a data frame containing at least the specified signal and relevant
%                 epocs for each trial.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       desired signal to be analyzed and epocs to mark each
%                       trial.
%
%       STREAMFIELDNAME:    The name of the field in the data structure that
%                       contains the signal to be cut.
%
%       startepocfieldname:      The name of the epoc to provide the start index for
%                       each trial.
%
%       endepocfieldname:        The name of the epoc to provide the end index for
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
    

% Find mean trial trace
    disp(['FINDING TRIAL MEAN: ', trialstreamfieldname])
    for eachfile = 1:length(data)
        trialdata = [data(eachfile).(trialstreamfieldname)];
        data(eachfile).(append(trialstreamfieldname,'_mean')) = mean(trialdata,1);
    end
end