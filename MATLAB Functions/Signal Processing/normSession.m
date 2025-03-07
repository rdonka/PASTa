function [data] = normSession(data,whichstream)
% NORMSESSION    Normalizes the data stream to zscore based on the whole session.
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:           Data structure; Must contain at least the stream 
%                       specified to be normalized.
%
%       WHICHSTREAM:    String; The name of the field containing the stream 
%                       to be normalized.
% OUTPUTS:
%       DATA:           Data structure; The original data structure with 
%                       'data.WHICHSTREAMz_normsession' added.
%
% Written by R M Donka, August 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Normalize to whole session
disp(append('NORM SESSION: Normalizing ',whichstream,' to whole session mean and standard deviation.'))
disp(append('   Normalized data will be output to the field: ',whichstream, 'z_normsession'))
    for eachfile = 1:length(data)
        disp(append('     NORMALIZING: File ',num2str(eachfile)))
        data(eachfile).(append(whichstream, 'z_normsession')) = (data(eachfile).(whichstream)-mean(data(eachfile).(whichstream),"omitnan"))/std(data(eachfile).(whichstream),"omitnan");
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