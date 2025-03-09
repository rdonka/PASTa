function [data] = cropFPdata(data,whichcropstart,whichcropend,whichstreams,varargin)
% CROPFPDATA    Crop specified data streams and adjust epocs for fiber photometry data.
%
%   CROPFPDATA(DATA, WHICHCROPSTART, WHICHCROPEND, WHICHSTREAMS, 'WHICHEPOCS', {'epoc1. 'epoc2'})
%   trims data streams between WHICHCROPSTART and WHICHCROPEND indices and 
%   adjusts specified epocs by the cropped amount.
%
% REQUIRED INPUTS:
%       DATA            - Struct array containing at least the specified 
%                         input fields.
%
%       WHICHCROPSTART  - String; Field name for cropping start index.
%                         Example: 'sessionstart'
%
%       WHICHCROPEND    - String; Field name for cropping end index.
%                         Example: 'sessionend'
%
%       WHICHSTREAMS    - Cell array of strings; Names of data streams to crop.
%                         Example: {'sig', 'baq'}
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%       'whichepocs'    - Cell array of epoc field names to adjust for cropping.
%                         Example: {'injt', 'trialstart', 'trialend'}
%
% OUTPUTS:
%       DATA            - Modified data structure with cropped streams
%                         and adjusted epocs (if specified).
%
% EXAMPLE:
%       data = cropFPdata(data, 'sessionstart', 'sessionend', {'sig', 'baq'}, ...
%                         'whichepocs', {'injt', 'trialstart'});
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/
    
    %% Prepare Settings
    % Import required and optional inputs into a structure
    p = inputParser;
    addParameter(p, 'whichepocs', {}, @(x) iscell(x) && all(cellfun(@ischar, x)));
    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Main display and function inputs
    if isempty(params.whichepocs)
        disp('CROPFPDATA: No epocs specified for adjustment.');
    else
        disp('CROPFPDATA: Specified epocs will be adjusted.');
        disp('   EPOCS:');
        disp(params.whichepocs);
    end

    % %% Prepare Settings
    % % Import required and optional inputs into a structure
    % inputs = struct(...
    %     'whichepocs',[]);
    % 
    % inputs = parseArgsLite(varargin,inputs);
    % 
    %     % Import required and optional inputs into a structure
    % p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    % addParameter(p, 'trim', defaultparameters.trim, @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'})); % trim: input must be a numeric value, nonnegative, and integer
    % addParameter(p, 'skipexisting', defaultparameters.skipexisting, @(x) any(x == [0 1])); % skipexisting: input must be 0 or 1
    % parse(p, varargin{:});

   

    %% Crop data
    for eachfile = 1:length(data) 
        disp(append('   Cropping File Number: ',num2str(eachfile)))
    
        startloc = data(eachfile).(whichcropstart); % Start index of the session
        endloc = data(eachfile).(whichcropend); % End index of the session
    
        % Loop through each specified stream and crop it
        for eachstream = 1:length(whichstreams)
            stream = char(whichstreams(eachstream)); % Get stream name as a string
            try 
                data(eachfile).(stream) = data(eachfile).(stream)(startloc:endloc);  % Crop the data based on start and end indices
            catch
                warning(['File ',num2str(eachfile),' - failed to crop stream: ',stream])  % Display warning if the stream cannot be cropped
            end
        end
        
        % Adjust specified epocs, if provided
        if isempty(params.whichepocs) == false
            whichepocs = params.whichepocs;
            for eachepoc = 1:length(whichepocs)
                epoc = char(whichepocs(eachepoc)); % Get epoc name
                try 
                    data(eachfile).(epoc) = data(eachfile).(epoc)-(startloc-1); % Adjust the epoc timing by shifting it based on the cropped start index
                catch
                    warning(['File ',num2str(eachfile),' - failed to adjust epoc: ',epoc]) % Display warning if the epoc cannot be adjusted
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