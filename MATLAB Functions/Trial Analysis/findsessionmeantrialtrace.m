function [data] = findsessionmeantrialtrace(data,trialstreamfieldname,trialtypefieldname,trialeventepocfieldnames,varargin)
% FINDSESSIONMEANTRIALTRACE    Finds the mean trace by trial type for each session. Requires individual 
%                              trial streams from cutTrialdata function and relevant event epochs. 
%                              Outputs mean trace for each session to an added field in the data structure.
%
% REQUIRED INPUTS:
%   DATA                        - Struct array; Must contain at least the desired cut trial stream 
%                                 to be normalized and epocs to mark the start and end of each trial's
%                                 baseline period.
%
%   TRIALSTREAMFIELDNAME        - String; The name of the field in the data structure that
%                                 contains the signal to be cut.
%
%   TRIALTYPEFIELDNAME          - String; The name of the field in the data structure containing trial
%                                 type IDs by which to average the data across. Must be spatially matched 
%                                 to the trial data streams in TRIALSTREAMFIELDNAME.
%
%   TRIALEVENTEPOCFIELDNAMES    - Cell array of strings; The names of the fields containing event epochs 
%                                 relevant for the trial analysis.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   TRIALIDOUPUTFIELDNAME       - String; Output name for the trial ID field output with the mean trace data. 
%                                 Default: 'TrialType'
% OUTPUT:
%   DATA    - Struct array; This is the original data structure with an added
%             field for the mean trial streams by trial type.
%               - data.mean_<trialstreamfieldname>: contains the mean trial data stream by trial type 
%                 and trial type IDs spatially matched to the stream data.
%                   - data.mean_<trialstreamfieldname>.<trialidoutputfieldname>:
%                     trial IDs
%                   - data.mean_<trialstreamfieldname>.trialdata: mean
%                   trial stream by trial type.
%                   -
%                   data.mean_<trialstreamfieldname>.<eventepochfieldname>:
%                   event epochs matched to trial type for field included
%                   in 'trialeventepocfieldnames'.
%
% Author:  Rachel Donka (2026)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

    %% Prepare Settings
    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'trialidoutputfieldname', 'TrialType'); 

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Initialize params
    params.trialstreamfieldname = trialstreamfieldname;
    params.trialtypefieldname = trialtypefieldname;
    params.trialeventepocfieldnames = trialeventepocfieldnames;
    
    % Prepare trial ID output field name
    if ~isempty(params.trialidoutputfieldname)
        trialidoutputfieldname = params.trialidoutputfieldname;
    end
    % Prepare trial stream name output field name
    outputtrialstreamfieldname = append('mean_',trialstreamfieldname);

    % Find trial stream means by trial type
    for eachfile = 1:length(data)
        fprintf('Find trial means for file number: %.f \n',eachfile) % Display which file is being subtracted
        trialtypes = unique([data(eachfile).(trialtypefieldname)]); % Identify unique trial types
        trialmeanstruct = struct(); % Prepare temp structure
        try
            % Find mean for each trial tupe
            for eachtrialtype = 1:length(trialtypes) % Save each trial type to a row
                currtrialtype = trialtypes(eachtrialtype);
                currtrialtypeidxs = find([data(eachfile).(trialtypefieldname)] == currtrialtype);

                trialmeanstruct.(trialidoutputfieldname)(eachtrialtype,1) = trialtypes(eachtrialtype); % Add trial type ID
                trialmeanstruct.trialdata(eachtrialtype,:) = mean(data(eachfile).(trialstreamfieldname)(currtrialtypeidxs,:)); % Add mean trial stream trace
                
                % Add event epochs
                for eachtrialeventepoc = 1:length(trialeventepocfieldnames)
                    currepocname = trialeventepocfieldnames{eachtrialeventepoc};
                    trialmeanstruct.(currepocname)(eachtrialtype,1) = mean(data(eachfile).(currepocname)(currtrialtypeidxs,:));
                end
            end

            data(eachfile).(outputtrialstreamfieldname) = trialmeanstruct; % Add subject trial means to data structure
        catch
            disp(['Find trial means without success:' num2str(eachfile)]);
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