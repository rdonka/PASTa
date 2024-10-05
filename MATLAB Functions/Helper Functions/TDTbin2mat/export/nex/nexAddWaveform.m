function [ nexFile ] = nexAddWaveform( nexFile, WFreq, timestamps, waveforms, name, preThresholdTimeInSeconds, numberOfPointsInWaveform, wireNumber, unitNumber )
% [nexFile] = nexAddWaveform( nexFile, WFreq, timestamps, waveforms, name, preThresholdTimeInSeconds, numberOfPointsInWaveform )
%             -- adds waveform variable to nexFile data structure
%
% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   startTime - time of the first data point in seconds
%   WFreq - A/D sampling rate of waveform variable in samples per second
%   timestamps - vector of wave timestamps (in seconds)
%   waveforms - matrix of waveform variable values in milliVolts
%               each waveform is a column (matrix can be empty if numberOfPointsInWaveform is specified)
%   name - waveform variable name  
%   preThresholdTimeInSeconds - (optional parameter) if waveform timestamp in seconds is t, 
%         then the timestamp of the first point of waveform is t - prethresholdTimeInSeconds
%   numberOfPointsInWaveform - (optional parameter) number of points in waveform; 
%         if this parameter is omitted, the number of points in waveform is 
%         the number of rows in waveforms matrix
%   wireNumber - wire number (optional parameter)
%   unitNumber - unit number (optional parameter). Zero unit number means unsorted.

    if WFreq > nexFile.freq
        error 'waveform frequency cannot exceed file timestamp frequency'
    end
    if isa(name, 'char') == 0
        error 'name should be of type char'
    end

    wire = 0;
    unit = 0;
    if nargin > 7
        wire = wireNumber;
    end
    if nargin > 8
        unit = unitNumber;
    end
    
    waveCount = 0;
    if(isfield(nexFile, 'waves'))
        waveCount = size(nexFile.waves, 1);
    end
    
    waveCount = waveCount+1;
    nexFile.waves{waveCount,1}.name = name;
    nexFile.waves{waveCount,1}.WFrequency = WFreq;
    nexFile.waves{waveCount,1}.wireNumber = wire;
    nexFile.waves{waveCount,1}.unitNumber = unit;
    
    if nargin < 6
        nexFile.waves{waveCount,1}.PrethresholdTimeInSeconds = 0;
        nexFile.waves{waveCount,1}.varVersion = 100;
    else
        nexFile.waves{waveCount,1}.PrethresholdTimeInSeconds = preThresholdTimeInSeconds;
        % set var version to 101 and file version to 106 to allow saving pre-threshold time
        % in .nex files
        nexFile.waves{waveCount,1}.varVersion = 102;
        nexFile.version = 106;
    end    
    
    if nargin < 7
        if size(waveforms, 1) == 0
           error 'if waveforms matrix is empty, specify number of points in waveform as the 7th parameter to nexAddWaveform' 
        end
        nexFile.waves{waveCount,1}.NPointsWave = size(waveforms, 1);
    else
        if numberOfPointsInWaveform <= 0
            error 'number of points in waveform should be positive'
        end
        nexFile.waves{waveCount,1}.NPointsWave = numberOfPointsInWaveform;
    end
    
    if size(waveforms,1) > 0 && size(waveforms,1) ~= nexFile.waves{waveCount,1}.NPointsWave
        error 'number of points in waveform is incorrect (should be equal to the number of rows in waveforms parameter)'
    end
    
    % timestamps should be a vector
    % if row, transpose a vector
    numTimestamps = size(timestamps,1);
    if size(timestamps,1) == 1
        numTimestamps = size(timestamps,2);
        % if row, transpose to vector
        nexFile.waves{waveCount,1}.timestamps = timestamps';
    else
        nexFile.waves{waveCount,1}.timestamps = timestamps;
    end
    if numTimestamps == size(waveforms, 2) 
        nexFile.waves{waveCount,1}.waveforms = waveforms;
    else
        error 'sizes of timestamps and waveforms do not match'
    end
    
    % modify end of file timestamp value in file header
    if size(nexFile.waves{waveCount,1}.timestamps,1) > 0 
        nexFile.tend = max(nexFile.tend, timestamps(end));
    end
end

