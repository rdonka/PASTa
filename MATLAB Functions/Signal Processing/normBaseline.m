function [data] = normBaseline(data,streamfieldname,BLstartfieldname,BLendfieldname)
% NORMBASELINE    Normalizes a specified data stream to z-score based on the 
%                 session baseline period.
%
%   NORMBASELINE(DATA, STREAMFIELDNAME, BLSTARTFIELDNAME, BLENDFIELDNAME) normalizes 
%   the data stream specified by STREAMFIELDNAME within the DATA structure to 
%   its z-score, using the mean and standard deviation calculated over the 
%   session baseline defined by BLSTARTFIELDNAME and BLENDFIELDNAME. The normalized 
%   data is added to the DATA structure with the field name 
%   '<streamfieldname>_z_normsession'.
%
% REQUIRED INPUTS:
%       DATA              - Structure array; each element represents a session and must
%                           contain the field specified by STREAMFIELDNAME.
%
%       STREAMFIELDNAME   - String; the name of the field within DATA to be normalized.
%
%       BLSTARTFIELDNAME  - String; The name of the field containing the index 
%                           of the start of the baseline period.
%
%       BLENDFIELDNAME    - String; The name of the field containing the index 
%                           of the end of the baseline period.
%
%   OUTPUTS:
%       DATA:           Structure array; the original DATA structure with an added
%                       field '<whichStream>z_normbaseline' containing the normalized data.
%
%   EXAMPLE:
%       % Assuming 'data' is a structure array:
%       data = normBaseline(data, 'sigfilt', 'blstart', 'blend');
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa

%% Normalize to session baseline
disp(['NORMBASELINE: Normalizing ',streamfieldname,' to session baseline mean and standard deviation.'])
disp(['   Baseline defined by ',BLstartfieldname,' and ',BLendfieldname,'.'])
disp(['   Normalized data will be output to the field: ',streamfieldname, 'z_normbaseline'])

    for eachfile = 1:length(data)
        disp(['   NORMALIZING: File ',num2str(eachfile)])
        try
            BLmean = mean(data(eachfile).(streamfieldname)(data(eachfile).(BLstartfieldname):data(eachfile).(BLendfieldname)),"omitnan"); % Find the mean of the session baseline
            BLsd = std(data(eachfile).(streamfieldname)(data(eachfile).(BLstartfieldname):data(eachfile).(BLendfieldname)),"omitnan"); % Find the standard deviation of the session baseline
    
            data(eachfile).(append(streamfieldname,'z_normbaseline')) = (data(eachfile).(streamfieldname) - BLmean)/BLsd; % Z score the whole session to the baseline
        catch
            warning(['File ',num2str(eachfile), ' - failed to normalize stream: ', streamfieldname]) 
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