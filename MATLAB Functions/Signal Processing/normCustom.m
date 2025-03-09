function [data] = normCustom(data,whichfullstream,whichcustomstream)
% NORMCUSTOM    Normalizes the whole session data stream based on a
%                custom period input as a separate stream field.
%
%   NORMBASELINE(DATA, WHICHSTREAM, WHICHBLSTART, WHICHBLEND) normalizes 
%   the data stream specified by WHICHSTREAM within the DATA structure to 
%   its z-score, using the mean and standard deviation calculated over the 
%   session baseline defined by WHICHBLSTART and WHICHBLEND. The normalized 
%   data is added to the DATA structure with the field name 
%   '<whichStream>_z_normsession'.
%
% REQUIRED INPUTS:
%   DATA                - Structure array; each element represents a session
%                         and must contain the fields specified by WHICHFULLSTREAM
%                         and WHICHCUSTOMSTREAM.
%
%   WHICHFULLSTREAM     - String; The name of the field containing the full 
%                         data stream to be normalized.
%
%   WHICHCUSTOMSTREAM   - String; The name of the field containing the cut 
%                         data stream to use as the reference for
%                         normalization.
%
%   OUTPUTS:
%       DATA:           Structure array; the original DATA structure with an added
%                       field '<whichStream>_z_normcustom' containing the normalized data.
%
%   EXAMPLE:
%       % Assuming 'data' is a structure array with a field 'sigfilt':
%       data = normCustom(data, 'sigfilt', 'customsigfilt');
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa


% INPUTS:
%       DATA:           Data structure; Must contain at least the stream to 
%                       be normalized, baseline start field, and baseline 
%                       end field.
%
%       WHICHFULLSTREAM: String; The name of the field containing the full 
%                       data stream to be normalized.
%
%       WHICHCUSTOMSTREAM: String; The name of the field containing the
%                       pre-prepared customized cut of the stream to be
%                       used as the reference for the normalization of the
%                       full stream. Example: 3 minutes before and after
%                       the trial portion of a session.
%
% OUTPUTS:
%       DATA:           Data structure; This is the original data structure 
%                       with 'data.WHICHSTREAMz_normcustom' added.
%
% Written by R M Donka, February 2025.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Normalize to session baseline
disp(append('NORM CUSTOM: Normalizing ',whichfullstream,' to mean and standard deviation of customized period of session.'))
disp(append('   Mean and SD defined by ',whichcustomstream,'.'))
disp(append('   Normalized data will be output to the field: ',whichfullstream, 'z_normcustom'))

    for eachfile = 1:length(data)
        disp(append('   NORMALIZING: File ',num2str(eachfile)))
        try
            BLmean = mean(data(eachfile).(whichcustomstream)); % Find the mean of the session baseline
            BLsd = std(data(eachfile).(whichcustomstream)); % Find the standard deviation of the session baseline
    
            data(eachfile).(append(whichfullstream,'z_normcustom')) = (data(eachfile).(whichfullstream) - BLmean)/BLsd; % Z score the whole session to the baseline
        catch
            disp(append('WARNING: File ',num2str(eachfile), ' - failed to normalize stream: ', whichfullstream)) 
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