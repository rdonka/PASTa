function [streamFFT,streamF] = preparestreamFFT(streamdata,fs)
% PREPARESTREAMFFT  Prepares the one sided FFT of the stream for plotting 
%                   and analysis.
% INPUTS:
%   STREAMDATA:         An array containing the values from a collected data
%                       stream.
% 
%   FS:                 The sampling rate of the stream in Hz. E.g., 1017.
%
% OPTIONAL INPUTS:
%
% OUTPUTS:
%       streamFFT:      An array containing the preapred FFT of the input 
%                       stream.
%
% Written by R M Donka, August 2024.
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.                

    L = length(streamdata);

% Take FFTs of signal and background
    fft1 = abs((fft(streamdata,L))/L); % Take the fft
    fft2 = fft1(1:floor(L/2)+1); % Take one side
    fft2(2:end-1) = 2*fft2(2:end-1); % Double the power to adjust for one sided
    streamFFT = fft2;

    streamF = fs*(0:floor(L/2))/L; % Spatially matched array frequency labels
end