function [data] = notchStreamArtifacts(data,streamfieldnames,notchstartfieldname,notchendfieldname)
% NOTCHSTREAMARTIFACTS  Removes artifacts from data streams using manually generated 
%                       notch start and notch end indexes.
%
%   DATA = NOTCHSTREAMARTIFACTS(DATA, STREAMFIELDNAMES, NOTCHSTARTFIELDNAME, NOTCHENDFIELDNAME)
%   detects artifacts based on amplitude deviations and replaces them with NaNs or mean values.
%
% REQUIRED INPUTS:
%       DATA              - Structure array; each element represents a session and must
%                           contain the fields specified by STREAMFIELDNAMES, 
%                           NOTCHSTARTFIELDNAME, and NOTCHENDFIELDNAME.
%
%       STREAMFIELDNAMES  - Cell array; Names of all streams for artifact
%                           removal.
%
%       NOTCHSTARTFIELDNAME - String; Name of field containing indexes for start 
%                             of artifact notch.
%
%       NOTCHENDFIELDNAME   - String; Name of field containing indexes for end 
%                             of artifact notch.
%
% OUTPUTS:
%       DATA - Struct with artifact notches replaced by NaNs for all
%       specified data streams.
%
%   EXAMPLE:
%       streamfieldnames = {'sig','baq','baqscaled','sigsub','sigfilt'}
%       notchstartfieldname = 'artifactstartidx'
%       notchendfieldname = 'artifactendidx'
%
%       [data] = notchstreamartifacts(data, streamfieldnames, notchstartfieldname, notchendfieldname);
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

% Display function inputs
disp('NOTCHSTREAMARTIFACTS: Stream artifacts will be removed for specified streams based on manually identified artifact start and end indexes.')
disp(['   STREAMFIELDNAMES: ', streamfieldnames])
disp(['   NOTCHSTARTFIELDNAME: ', notchstartfieldname])
disp(['   NOTCHENDFIELDNAME: ', notchendfieldname])  


    %% Remove Artifacts
    for eachfile = 1:length(data)
        if isempty(data(eachfile).(notchstartfieldname))
            continue
        else
            for eachstream = 1:length(streamfieldnames)
                currstreamfieldname = streamfieldnames{eachstream};
                for eachnotch = 1:length(data(eachfile).(notchstartfieldname))
                    currnotchstart = data(eachfile).(notchstartfieldname)(eachnotch);
                    currnotchend = data(eachfile).(notchendfieldname)(eachnotch);
                    data(eachfile).(currstreamfieldname)(currnotchstart:currnotchend) = NaN;
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
