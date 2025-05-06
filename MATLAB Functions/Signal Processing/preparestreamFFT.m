function [streamFFT, streamF] = preparestreamFFT(streamdata,fs)
% PREPARESTREAMFFT  Compute the one-sided FFT of a data stream for analysis and plotting.
%
%   [streamFFT, streamF] = PREPARESTREAMFFT(streamdata, fs) computes the
%   one-sided Fast Fourier Transform (FFT) of the input data stream and
%   returns the FFT magnitudes along with their corresponding frequencies.
%
% REQUIRED INPUTS:
%       STREAMDATA - Numeric vector; the input data stream to be transformed.
%
%       FS         - Positive scalar; the sampling rate of the data stream
%                    in Hz.
%
%   OUTPUTS:
%       STREAMFFT  - Numeric vector; the magnitudes of the one-sided FFT
%                    of the input data stream.
%
%       STREAMF    - Numeric vector; the corresponding frequency values in
%                    Hz for the FFT magnitudes.
%
%   EXAMPLE:
%       % Compute the FFT
%       [streamFFT, streamF] = preparestreamFFT(streamdata, fs);
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

    %% Compute FFT
    L = length(streamdata); % Length of data
    fft1 = abs((fft(streamdata,L))/L); % Take the fft
    fft2 = fft1(1:floor(L/2)+1); % Take one side
    fft2(2:end-1) = 2*fft2(2:end-1); % Double the power to adjust for one-sided
    
    streamFFT = fft2; % Output streamFFT
    streamF = fs*(0:floor(L/2))/L; % Output streamF; Spatially matched array frequency labels
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