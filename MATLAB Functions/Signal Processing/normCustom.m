function [data] = normCustom(data,fullstreamfieldname,customstreamfieldname,varargin)
% NORMCUSTOM    Normalizes the whole session data stream based on a
%                custom period input as a separate stream field.
%
%   NORMCUSTOM(DATA, FULLSTREAMFIELDNAME, CUSTOMSTREAMFIELDNAME) normalizes 
%   the data stream specified by FULLSTREAMFIELDNAME within the DATA structure 
%   to its z-score, using the mean and standard deviation calculated over the 
%   specified custom stream cut specified by CUSTOMSTREAMFIELDNAME in the 
%   DATA structure. The normalized data is added to the DATA structure with 
%   the field name '<fullstreamfieldname>_z_normcustom'.
%
% REQUIRED INPUTS:
%   DATA                   - Structure array; each element represents a session
%                            and must contain the fields specified by
%                            FULLSTREAMFIELDNAME and CUSTOMSTREAMFIELDNAME.
%
%   FULLSTREAMFIELDNAME   - String; The name of the field containing the full 
%                           data stream to be normalized.
%
%   CUSTOMSTREAMFIELDNAME - String; The name of the field containing the cut 
%                           data stream to use as the reference for
%                           normalization.
%
% OPTIONAL INPUTS:
%   NORMFIELDOUTPUTNAME:  - String; output name to be appended to the FULLSTREAMFIELDNAME 
%                           for the normalized stream. 
%                           Default: <fullstreamfieldname>_normcustom.
% OUTPUTS:
%   DATA:   Structure array; the original DATA structure with an added field 
%           '<fullstreamfieldname>z_normcustom' containing the normalized data.
%
%   EXAMPLE:
%       % Assuming 'data' is a structure array with a field 'sigfilt':
%       data = normCustom(data, 'sigfilt', 'customsigfilt');
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa

% Prepare default values
defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.

% Import required and optional inputs into a structure
p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details

% Add optional name-value pair arguments with validation
addParameter(p, 'normfieldoutputname', defaultparameters.normfieldoutputname, @(x) ischar(x)); % normfieldoutputname: input must be char

parse(p, varargin{:});

% Retrieve all parsed inputs into params structure
allparams = p.Results;

normfieldoutputname = append(fullstreamfieldname,'z_',allparams.normfieldoutputname);

%% Normalize to session baseline
disp(['NORM CUSTOM: Normalizing ',fullstreamfieldname,' to mean and standard deviation of customized period of session.'])
disp(['   Mean and SD defined by ',customstreamfieldname,'.'])
disp(['   Normalized data will be output to the field: ',normfieldoutputname])

    for eachfile = 1:length(data)
        disp(['   NORMALIZING: File ',num2str(eachfile)])
        try
            BLmean = mean(data(eachfile).(customstreamfieldname),'omitnan'); % Find the mean of the session baseline
            BLsd = std(data(eachfile).(customstreamfieldname),'omitnan'); % Find the standard deviation of the session baseline
    
            data(eachfile).(normfieldoutputname) = (data(eachfile).(fullstreamfieldname) - BLmean)/BLsd; % Z score the whole session to the baseline
            data(eachfile).(append(normfieldoutputname,'_mean')) = BLmean; % Add mean used for Z score to data structure
            data(eachfile).(append(normfieldoutputname,'_sd')) = BLsd; % Add sd used for Z score to data structure
        catch
            warning(['File ',num2str(eachfile), ' - failed to normalize stream: ', fullstreamfieldname]) 
        end
    end
end

% Copyright (C) 2024 Rachel Donka
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