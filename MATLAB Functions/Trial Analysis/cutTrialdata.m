function [data] = cutTrialdata(data,streamfieldname,startepocfieldname,endepocfieldname,varargin)
% CUTTRIALDATA    Cuts fiber photometry signal data into individual trials based on epocs.
%                 Requires a data structure containing at least the specified signal and relevant
%                 epocs for each trial.
%
% REQUIRED INPUTS:
%       DATA              - Structure array; each element represents a session and must
%                           contain the fields specified by STREAMFIELDNAME, STARTEPOCFIELDNAME, 
%                           and ENDEPOCFIELDNAME.
%
%       STREAMFIELDNAME   - String; the name of the field within DATA to be normalized.
%
%       STARTEPOCFIELDNAME  - String; The name of the field containing the index 
%                           of the start of the trials.
%
%       ENDEPOCFIELDNAME    - String; The name of the field containing the index 
%                           of the end of the trials.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'epocsfieldnames'    - Cell array of strings; names of additional event epochs to be adjusted 
%                          for temporal alignment with cut trial data streams.
%   OUTPUTS:
%       DATA:   - Structure array; the original DATA structure with added fields:
%                   - 'trial_<whichStream>' containing the cut trial data streams.
%                   - 'trial_<STARTEPOCFIELDNAME>' containing the adjusted
%                      trial start index.
%                   - 'trial_<STARTEPOCFIELDNAME>' containing the adjusted
%                      trial end index.
%                   - 'trial_<EPOCSFIELDNAMES>' fields for each additional event epoch
%                      adjusted to the indexes that correspond to the cut
%                      trial data streams.
%
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
% EXAMPLE USAGE:
%   [data] = cuttrialdata(data, 'sigfilt','trialstart', 'trialend', ... % Required inputs
%                           'epocsfieldnames', {'infonset','infoffset'}) % Optional additional epocs to adjust
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'epocsfieldnames', ''); % exportfilepath: string with the custom exportfilepath

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Cut trial data
        disp(append('CUTTING TRIAL DATA: ', streamfieldname, ' by indices in ', startepocfieldname, ' and ', endepocfieldname))
        for eachfile = 1:length(data)
            maxtrialsamples = max(data(eachfile).(endepocfieldname) - data(eachfile).(startepocfieldname)) + 1;
            ntrials = length(data(eachfile).(startepocfieldname));
            % Prepare data structure
            data(eachfile).(append('trial_',streamfieldname))(1:ntrials,1:maxtrialsamples) = NaN;
    
            for eachtrial = 1:size(data(eachfile).(startepocfieldname),1)
                currtrialsamples = data(eachfile).(endepocfieldname)(eachtrial) - data(eachfile).(startepocfieldname)(eachtrial)+1;
                data(eachfile).(append('trial_',streamfieldname))(eachtrial,1:currtrialsamples) = data(eachfile).(streamfieldname)(data(eachfile).(startepocfieldname)(eachtrial):data(eachfile).(endepocfieldname)(eachtrial));
                data(eachfile).(append('trial_',startepocfieldname))(eachtrial,1) = 1;
                data(eachfile).(append('trial_',endepocfieldname))(eachtrial,1) = data(eachfile).(endepocfieldname)(eachtrial) - data(eachfile).(startepocfieldname)(eachtrial);
                if isempty(params.epocsfieldnames) % Adjust additional inputs if specified with optional input WHICHEPOCS
                    continue
                else
                    for eachepoc = 1:length(params.epocsfieldnames)
                        epoc = char(params.epocsfieldnames(eachepoc));
                        data(eachfile).(append('trial_',epoc))(eachtrial,1) = data(eachfile).(epoc)(eachtrial)-data(eachfile).(startepocfieldname)(eachtrial)+1;
                    end
                end
            end
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