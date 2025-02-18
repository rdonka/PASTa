function [data] = centerTrialdata(data,whichsignal,startepoc,endepoc,varargin)
% CUTTRIALDATA    Cuts fiber photometry signal data into individual trials based on epocs.
%                 Requires a data frame containing at least the specified signal and relevant
%                 epocs for each trial.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       desired signal to be analyzed and epocs to mark each
%                       trial.
%
%       WHICHSIGNAL:    The name of the field in the data structure that
%                       contains the signal to be cut.
%
%       STARTEPOC:      The name of the epoc to provide the start index for
%                       each trial.
%
%       ENDEPOC:        The name of the epoc to provide the end index for
%                       each trial.
%
%       WHICHEPOCS:      A list of the names of the fields of any other epocs
%                       relevant to trial analysis.
% OUTPUTS:
%       DATA:           This is the original data structure with added
%                       fields for trial signal data and epocs.
%
%                       data.trial_WHICHSIGNAL: contains the cut trial data
%                       with each trial saved as a row.
%
%                       data.trial_startepoc: contains conthe index for the epoc
%                       specified as epoc start relative to the trimmed
%                       trial data.
%
%                       data.trial_whichepocs: each other epoc will generate
%                       a field that contains the index for that epoc
%                       relative to the startepoc for each trial.
%
% Written by R M Donka, January 2024
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

% Set up optional inputs
optionalinputs = struct('whichepocs',[]);
optionalinputs = parseArgsLite(varargin,optionalinputs);
    

% Cut trial data
    disp(append('CUTTING TRIAL DATA: ', whichsignal, ' by indices in ', startepoc, ' and ', endepoc))
    for eachfile = 1:length(data)
        for eachtrial = 1:size(data(eachfile).(startepoc),1)
            data(eachfile).(append('trial_',whichsignal))(eachtrial,:) = data(eachfile).(whichsignal)(data(eachfile).(startepoc)(eachtrial):data(eachfile).(endepoc)(eachtrial));
            data(eachfile).(append('trial_',startepoc))(eachtrial,1) = 1;
            data(eachfile).(append('trial_',endepoc))(eachtrial,1) = data(eachfile).(endepoc)(eachtrial) - data(eachfile).(startepoc)(eachtrial);
            if isempty(optionalinputs.whichepocs) % Adjust additional inputs if specified with optional input WHICHEPOCS
                continue
            else
                for eachepoc = 1:length(optionalinputs.whichepocs)
                    epoc = char(optionalinputs.whichepocs(eachepoc));
                    data(eachfile).(append('trial_',epoc))(eachtrial,1) = data(eachfile).(epoc)(eachtrial)-data(eachfile).(startepoc)(eachtrial);
                end
            end
        end
    end
end