function [streamFFT,streamF] = preparestreamFFT(streamdata,fs)
% PREPARESTREAMFFT  Prepares the one sided FFT of the stream for plotting 
%                   and analysis.
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%   STREAMDATA:         Array; Values from a single collected data stream.
%
%   FS:                 Sampling rate of the data stream.
%
% OUTPUTS:
%   STREAMFFT:          An array containing the prepared FFT magnitudes of
%                       the input stream.
%
%   STREAMF:            Corresponding frequencies to the prepared FFT
%                       magnitudes in STREAMFFT.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

    %% Find length of data
    L = length(streamdata);

    %% Prepare the FFT
    fft1 = abs((fft(streamdata,L))/L); % Take the fft
    fft2 = fft1(1:floor(L/2)+1); % Take one side
    fft2(2:end-1) = 2*fft2(2:end-1); % Double the power to adjust for one-sided
    streamFFT = fft2; % Output streamFFT
    streamF = fs*(0:floor(L/2))/L; % Output streamF; Spatially matched array frequency labels
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