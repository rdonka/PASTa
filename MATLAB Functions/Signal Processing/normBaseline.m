function [data] = normBaseline(data,whichstream,whichblstart,whichblend)
% NORMSESSION    Normalizes the whole session data stream based on a
%                baseline period.
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:           Data structure; Must contain at least the stream to 
%                       be normalized, baseline start field, and baseline 
%                       end field.
%
%       WHICHSTREAM:    String; The name of the field containing the data 
%                       stream to be normalized.
%
%       WHICHBLSTART:   String; The name of the field containing the index 
%                       of the start of the baseline period.
%
%       WHICHBLEND:     String; The name of the field containing the index 
%                       of the end of the baseline period.
%
% OUTPUTS:
%       DATA:           Data structure; This is the original data structure 
%                       with 'data.WHICHSTREAMz_normbaseline' added.
%
% Written by R M Donka, August 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Normalize to session baseline
disp(append('NORM BASELINE: Normalizing ',whichstream,' to session baseline mean and standard deviation.'))
disp(append('     Baseline defined by ',whichblstart,' and ',whichblend,'.'))

    for eachfile = 1:length(data)
        disp(append('   Normalizing: File ',num2str(eachfile)))
        try
            BLmean = mean(data(eachfile).(whichstream)(data(eachfile).(whichblstart):data(eachfile).(whichblend))); % Find the mean of the session baseline
            BLsd = std(data(eachfile).(whichstream)(data(eachfile).(whichblstart):data(eachfile).(whichblend))); % Find the standard deviation of the session baseline
    
            data(eachfile).(append(whichstream,'z_normbaseline')) = (data(eachfile).(whichstream) - BLmean)/BLsd; % Z score the whole session to the baseline
        catch
            disp(append('WARNING: File ',num2str(eachfile), ' - failed to normalize stream: ', whichstream)) 
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