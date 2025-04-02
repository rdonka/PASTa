function [data] = normSession(data,whichstream)
% NORMSESSION    Normalizes a specified data stream to z-score based on the entire session.
%
%   NORMSESSION(DATA, WHICHSTREAM) normalizes the data stream specified by WHICHSTREAM
%   within the DATA structure to its z-score, using the mean and standard deviation
%   calculated over the entire session. The normalized data is added to the DATA
%   structure with the field name '<whichStream>_z_normsession'.
%
% REQUIRED INPUTS:
%       DATA          - Structure array; each element represents a session
%                       and must contain the field specified by WHICHSTREAM.
%
%       WHICHSTREAM   - String; the name of the field within DATA to be normalized.
%
%   OUTPUTS:
%       DATA           - Structure array; the original DATA structure with an added
%                       field '<whichStream>_z_normsession' containing the normalized data.
%
%   EXAMPLE:
%       % Assuming 'data' is a structure array with a field 'sigfilt':
%       data = normSession(data, 'sigfilt');
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Normalize to whole session mean and SD
disp(['NORMSESSION: Normalizing ',whichstream,' to whole session mean and standard deviation.'])
disp(['   Normalized data will be output to the field: ',whichstream, 'z_normsession'])
    for eachfile = 1:length(data)
        try
            disp(append('     NORMALIZING: File ',num2str(eachfile)))
            data(eachfile).(append(whichstream, 'z_normsession')) = (data(eachfile).(whichstream)-mean(data(eachfile).(whichstream),"omitnan"))/std(data(eachfile).(whichstream),"omitnan");
         catch
            warning(['File ',num2str(eachfile), ' - failed to normalize stream: ', whichstream]) 
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