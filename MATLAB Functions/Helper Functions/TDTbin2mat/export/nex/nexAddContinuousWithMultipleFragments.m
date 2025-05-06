function [ nexFile ] = nexAddContinuousWithMultipleFragments( nexFile, fragmentTimestamps, fragmentIndexes, adFreq, values, name )
% [nexFile] = nexAddContinuousWithMultipleFragments( nexFile, fragmentTimestamps, fragmentIndexes, adFreq, values, name )
%         -- adds continuous variable to nexFile data structure
%
% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   startTimes - array of fragment start times in seconds
%   fragmentStarts - array of fragment start times in seconds
%   adFreq - A/D sampling rate of continuous variable in samples per second
%   values - vector of continuous variable values in milliVolts
%   name - continuous variable name  
% 
    if adFreq > nexFile.freq
        error 'continuous sampling rate cannot exceed file timestamp frequency'
    end
    if isa(name, 'char') == 0
        error 'name should be of type char'
    end

    contCount = 0;
    if(isfield(nexFile, 'contvars'))
        contCount = size(nexFile.contvars, 1);
    end
    contCount = contCount+1;
    nexFile.contvars{contCount,1}.name = name;
    nexFile.contvars{contCount,1}.varVersion = 100;
    nexFile.contvars{contCount,1}.ADFrequency = adFreq;
    nexFile.contvars{contCount,1}.timestamps = fragmentTimestamps;
    nexFile.contvars{contCount,1}.fragmentStarts = fragmentIndexes;
    nexFile.contvars{contCount,1}.data = values;
    
    % values should be a vector
    if size(values,1) == 1
        % if row, transpose to vector
        nexFile.contvars{contCount,1}.data = values';
    else
        nexFile.contvars{contCount,1}.data = values;
    end
    
    % modify end of file timestamp value in file header
    % number of data points in the last fragment
    if size(nexFile.contvars{contCount,1}.timestamps,1) > 0 
        pointsInLastFragment = length(values) - fragmentIndexes(end) + 1;
        nexFile.tend = max(nexFile.tend, fragmentTimestamps(end)+(pointsInLastFragment - 1)/adFreq);
    end    
end
