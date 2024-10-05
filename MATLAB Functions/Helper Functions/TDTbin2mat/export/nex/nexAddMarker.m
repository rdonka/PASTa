function [ nexFile ] = nexAddMarker( nexFile, timestamps, fieldName, fieldValues, name )
% [nexFile] = nexAddMarker( nexFile, timestamps, name ) -- adds a marker event 
%             (series of timestamps with string tags) to nexFile data structure.
%             In general, each marker timestamp can have several tags.
%             For simplicity, this function allows to add a marker with a
%             single string tag for each timestamp.
%
% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   timestamps - vector of marker timestamps in seconds
%   fieldName - marker field name
%   fieldValues - marker field values as a cell array of strings, 
%      for example, {'011','012','013'}
%   name - marker name
    
    if numel(timestamps) ~= numel(fieldValues)
        error('different number of timestamps and field values');
    end
    if isa(name, 'char') == 0
        error 'name should be of type char'
    end

    % check that all field values are strings
    for i=1:numel(fieldValues)
        if ~ischar(fieldValues{i})
            error('field values should be strings');
        end
    end
    
    markerCount = 0;
    if(isfield(nexFile, 'markers'))
        markerCount = size(nexFile.markers, 1);
    end
    markerCount = markerCount+1;
    nexFile.markers{markerCount,1}.name = name;
    nexFile.markers{markerCount,1}.varVersion = 100;
    nexFile.markers{markerCount,1}.values{1,1}.name = fieldName;
    
    if size(fieldValues,1) == 1
        % if row, transpose to vector
        nexFile.markers{markerCount,1}.values{1,1}.strings = fieldValues';
    else
        nexFile.markers{markerCount,1}.values{1,1}.strings = fieldValues;
    end
    
    % timestamps should be a vector
    if size(timestamps,1) == 1
        % if row, transpose to vector
        nexFile.markers{markerCount,1}.timestamps = timestamps';
    else
        nexFile.markers{markerCount,1}.timestamps = timestamps;
    end

    % modify end of file timestamp value in file header
    if size(nexFile.markers{markerCount,1}.timestamps,1) > 0
        nexFile.tend = max(nexFile.tend, timestamps(end));
    end
end
