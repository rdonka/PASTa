function [streamFFT,streamF] = preparestreamFFT(streamdata,fs)
% PREPARESTREAMFFT  Prepares the one sided FFT of the stream for plotting 
%                   and analysis.
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
% Written by R M Donka, August 2024
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