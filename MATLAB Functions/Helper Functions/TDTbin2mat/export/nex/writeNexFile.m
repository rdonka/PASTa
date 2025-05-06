function [result] = writeNexFile(nexFile, fileName)
% [result] = writeNexFile(nexFile, fileName) -- write nexFile structure
% to the specified .nex file. returns 1 if succeeded, 0 if failed.
%
% INPUT:
%   nexFile - a structure containing .nex file data
%
%           SOME FIELDS OF THIS STRUCTURE (VERSIONS ETC.) ARE NOT DESCRIBED
%           BELOW. IT IS RECOMMENDED THAT YOU READ A VALID .NEX FILE
%           TO FILL THIS STRUCTURE, THEN MODIFY THE STRUCTURE AND SAVE IT.
%
%           IF YOU WANT TO CREATE NEW .NEX FILE, USE nexCreateFileData.m,
%           nexAddContinuous.m etc. See exampleSaveDataInNexFile.m.
%
%   fileName - if empty string, will use File Save dialog
%
%   nexFile - a structure containing .nex file data
%   nexFile.version - file version
%   nexFile.comment - file comment
%   nexFile.tbeg - beginning of recording session (in seconds)
%   nexFile.tend - end of recording session (in seconds)
%
%   nexFile.neurons - array of neuron structures
%           neurons{i}.name - name of a neuron variable
%           neurons{i}.timestamps - array of neuron timestamps (in seconds)
%               to access timestamps for neuron 2 use {n} notation:
%               nexFile.neurons{2}.timestamps
%
%   nexFile.events - array of event structures
%           events{i}.name - name of event variable
%           events{i}.timestamps - array of event timestamps (in seconds)
%
%   nexFile.intervals - array of interval structures
%           intervals{i}.name - name of interval variable
%           intervals{i}.intStarts - array of interval starts (in seconds)
%           intervals{i}.intEnds - array of interval ends (in seconds)
%
%   nexFile.waves - array of wave structures
%           waves{i}.name - name of waveform variable
%           waves{i}.NPointsWave - number of data points in each wave
%           waves{i}.WFrequency - A/D frequency for wave data points
%           waves{i}.timestamps - array of wave timestamps (in seconds)
%           waves{i}.waveforms - matrix of waveforms (in milliVolts), each
%                             waveform is a column
%
%   nexFile.contvars - array of continuous variable structures
%           contvars{i}.name - name of continuous variable
%           contvars{i}.ADFrequency - A/D frequency for data points
%
%           Continuous (a/d) data for one channel is allowed to have gaps
%           in the recording (for example, if recording was paused, etc.).
%           Therefore, continuous data is stored in fragments.
%           Each fragment has a timestamp and an index of the first data
%           point of the fragment (data values for all fragments are stored
%           in one array and the index indicates the start of the fragment
%           data in this array).
%           The timestamp corresponds to the time of recording of
%           the first a/d value in this fragment.
%
%           contvars{i}.timestamps - array of timestamps (fragments start times in seconds)
%           contvars{i}.fragmentStarts - array of start indexes for fragments in contvar.data array
%           contvars{i}.data - array of data points (in milliVolts)
%
%   nexFile.popvectors - array of population vector structures
%           popvectors{i}.name - name of population vector variable
%           popvectors{i}.weights - array of population vector weights
%
%   nexFile.markers - array of marker structures
%           markers{i}.name - name of marker variable
%           markers{i}.timestamps - array of marker timestamps (in seconds)
%           markers{i}.values - array of marker value structures
%               markers{i}.value.name - name of marker value
%               markers{i}.value.strings - array of marker value strings
%

result = 0;

if (nargin < 2 || isempty(fileName))
    [fname, pathname] = uiputfile('*.nex', 'Save as');
    if isequal(fname,0)
        error 'File name was not selected'
    end
    fileName = fullfile(pathname, fname);
end

% note 'l' option when opening the file.
% this options means that the file is 'little-endian'.
% this should ensure that the files are written correctly
% on big-endian systems, such as Mac G5.
fid = fopen(fileName, 'w', 'l', 'US-ASCII');
if(fid == -1)
    error 'Unable to open file'
end

maxTimestampTick = maxTimestampInTicks(nexFile);

if maxTimestampTick > pow2(31)
    error 'Unable to save .nex file: maximum timestamp exceeds 32-bit max. You may save data in .nex5 instead using writeNex5File function'
end

% count all the variables
neuronCount = 0;
eventCount = 0;
intervalCount = 0;
waveCount = 0;
contCount = 0;
markerCount = 0;

if(isfield(nexFile, 'neurons'))
    neuronCount = size(nexFile.neurons, 1);
end
if(isfield(nexFile, 'events'))
    eventCount = size(nexFile.events, 1);
end
if(isfield(nexFile, 'intervals'))
    intervalCount = size(nexFile.intervals, 1);
end
if(isfield(nexFile, 'waves'))
    waveCount = size(nexFile.waves, 1);
end
if(isfield(nexFile, 'contvars'))
    contCount = size(nexFile.contvars, 1);
end
if(isfield(nexFile, 'markers'))
    markerCount = size(nexFile.markers, 1);
end

nvar = int32(neuronCount+eventCount+intervalCount+waveCount+contCount+markerCount);

% write header information
fwrite(fid, 827868494, 'int32');
if(isfield(nexFile, 'version'))
    fwrite(fid, nexFile.version, 'int32');
else
    fwrite(fid, 106, 'int32');
end

writeStringPaddedWithZeros(fid, nexFile.comment, 256);
fwrite(fid, nexFile.freq, 'double');
fwrite(fid, int32(nexFile.tbeg*nexFile.freq), 'int32');
fwrite(fid, int32(maxTimestampTick), 'int32');
fwrite(fid, nvar, 'int32');

% skip location of next header and padding
fwrite(fid, char(zeros(1, 260)), 'char');

% calculate where variable data starts
dataOffset = int64(544 + nvar*208);

% write variable headers

varHeader.Type = 0;
varHeader.Version = 100;
varHeader.Name = ' ';
varHeader.DataOffset = 0;
varHeader.Count = 0;
varHeader.WireNumber = 0;
varHeader.UnitNumber = 0;
varHeader.Gain = 0;
varHeader.Filter = 0;
varHeader.XPos = 0;
varHeader.YPos = 0;
varHeader.WFrequency = 0;
varHeader.ADtoMV = 0;
varHeader.NPointsWave = 0;
varHeader.NMarkers = 0;
varHeader.MarkerLength = 0;
varHeader.MVOffset = 0;
varHeader.PrethresholdTimeInSeconds = 0;

% write neuron headers
for i = 1:neuronCount
    varHeader.Type = 0;
    varHeader.Version = 100;
    varHeader.Name = nexFile.neurons{i}.name;
    varHeader.Count = size(nexFile.neurons{i}.timestamps,1);
    varHeader.DataOffset = dataOffset;
    varHeader.WireNumber = 0;
    varHeader.UnitNumber = 0;
    varHeader.XPos = 0;
    varHeader.YPos = 0;
    
    if(isfield(nexFile.neurons{i}, 'varVersion'))
        varHeader.Version = nexFile.neurons{i}.varVersion;
    end
    
    if(isfield(nexFile.neurons{i}, 'wireNumber'))
        varHeader.WireNumber = nexFile.neurons{i}.wireNumber;
        varHeader.Version = max(101, varHeader.Version);
    end
    
    if(isfield(nexFile.neurons{i}, 'unitNumber'))
        varHeader.UnitNumber = nexFile.neurons{i}.unitNumber;
        varHeader.Version = max(101, varHeader.Version);
    end
    
    if(isfield(nexFile.neurons{i}, 'xPos'))
        varHeader.XPos = nexFile.neurons{i}.xPos;
    end
    if(isfield(nexFile.neurons{i}, 'yPos'))
        varHeader.YPos = nexFile.neurons{i}.yPos;
    end
    
    writeNexVarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*4;
    if dataOffset >= 2147483647
        fclose(fid);
        delete(fileName);
        error 'Data will not fit into .nex file (2 GB limit). Save data as .nex5 file instead.'
    end
end

% event headers
varHeader.Version = 100;
varHeader.WireNumber = 0;
varHeader.UnitNumber = 0;
varHeader.XPos = 0;
varHeader.YPos = 0;

for i = 1:eventCount
    varHeader.Type = 1;
    varHeader.Name = nexFile.events{i}.name;
    varHeader.Count = size(nexFile.events{i}.timestamps,1);
    varHeader.DataOffset = dataOffset;
    
    writeNexVarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*4;
    if dataOffset >= 2147483647
        fclose(fid);
        delete(fileName);
        error 'Data will not fit into .nex file (2 GB limit). Save data as .nex5 file instead.'
    end
end

% interval headers
for i = 1:intervalCount
    % interval variable type is 2
    varHeader.Type = 2;
    varHeader.Name = nexFile.intervals{i}.name;
    varHeader.Count = size(nexFile.intervals{i}.intStarts,1);
    varHeader.DataOffset = dataOffset;
    
    writeNexVarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*8;
    if dataOffset >= 2147483647
        fclose(fid);
        delete(fileName);
        error 'Data will not fit into .nex file (2 GB limit). Save data as .nex5 file instead.'
    end
end

% wave headers
for i = 1:waveCount
    varHeader.Count = size(nexFile.waves{i}.timestamps,1);
    c = 1;
    if varHeader.Count > 0
        % we need to recalculate a/d to millivolts factor
        wmin = min(min(nexFile.waves{i}.waveforms));
        wmax = max(max(nexFile.waves{i}.waveforms));
        c = max(abs(wmin),abs(wmax));
        if (c == 0)
            c = 1;
        else
            c = c/32767;
        end
    end
    nexFile.waves{i}.ADtoMV = c;
    nexFile.waves{i}.MVOffset = 0;
    
    % wave variable type is 3
    varHeader.Type = 3;
    varHeader.Version = nexFile.waves{i}.varVersion;
    varHeader.Name = nexFile.waves{i}.name;
    varHeader.DataOffset = dataOffset;
    
    varHeader.WireNumber = 0;
    varHeader.UnitNumber = 0;
    
    varHeader.WFrequency = nexFile.waves{i}.WFrequency;
    varHeader.ADtoMV = nexFile.waves{i}.ADtoMV;
    varHeader.NPointsWave = nexFile.waves{i}.NPointsWave;
    varHeader.MVOffset = nexFile.waves{i}.MVOffset;
    
    if(isfield(nexFile.waves{i}, 'varVersion'))
        varHeader.Version = nexFile.waves{i}.varVersion;
    end
    
    if(isfield(nexFile.waves{i}, 'wireNumber'))
        varHeader.WireNumber = nexFile.waves{i}.wireNumber;
        varHeader.Version = max(101, varHeader.Version);
    end
    
    if(isfield(nexFile.waves{i}, 'unitNumber'))
        varHeader.UnitNumber = nexFile.waves{i}.unitNumber;
        varHeader.Version = max(101, varHeader.Version);
    end
    
    if(isfield(nexFile.waves{i}, 'PrethresholdTimeInSeconds'))
        varHeader.PrethresholdTimeInSeconds = nexFile.waves{i}.PrethresholdTimeInSeconds;
        varHeader.Version = max(102, varHeader.Version);
    end
    
    writeNexVarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*4 + varHeader.NPointsWave*varHeader.Count*2;
    if dataOffset >= 2147483647
        fclose(fid);
        delete(fileName);
        error 'Data will not fit into .nex file (2 GB limit). Save data as .nex5 file instead.'
    end
end

% continuous variables
varHeader.Version = 100;
varHeader.WireNumber = 0;
varHeader.UnitNumber = 0;
varHeader.PrethresholdTimeInSeconds = 0;

for i = 1:contCount
    varHeader.Count = size(nexFile.contvars{i}.timestamps,1);
    
    c = 1;
    
    if varHeader.Count > 0
        wmin = min(min(nexFile.contvars{i}.data));
        wmax = max(max(nexFile.contvars{i}.data));
        c = max(abs(wmin),abs(wmax));
        if (c == 0)
            c = 1;
        else
            c = c/32767;
        end
    end
    nexFile.contvars{i}.ADtoMV = c;
    nexFile.contvars{i}.MVOffset = 0;
    
    % cont. variable type is 5
    varHeader.Type = 5;
    varHeader.Version = 100;
    varHeader.Name = nexFile.contvars{i}.name;
    
    varHeader.DataOffset = dataOffset;
    
    varHeader.WFrequency = nexFile.contvars{i}.ADFrequency;
    varHeader.ADtoMV = nexFile.contvars{i}.ADtoMV;
    varHeader.NPointsWave = size(nexFile.contvars{i}.data, 1);
    varHeader.MVOffset = nexFile.contvars{i}.MVOffset;
    
    writeNexVarHeader(fid, varHeader);
    % we have timestamps and indexes, so we use 8 bytes per count
    dataOffset = dataOffset + varHeader.Count*8 + varHeader.NPointsWave*2;
    if dataOffset >= 2147483647
        fclose(fid);
        delete(fileName);
        error 'Data will not fit into .nex file (2 GB limit). Save data as .nex5 file instead.'
    end
end

% markers
varHeader.WFrequency = 0;
varHeader.ADtoMV = 0;
varHeader.NPointsWave = 0;
varHeader.NMarkers = 0;
varHeader.MarkerLength = 0;
varHeader.MVOffset = 0;

for i = 1:markerCount
    
    nexFile.markers{i}.NMarkers = size(nexFile.markers{i}.values, 1);
    nexFile.markers{i}.MarkerLength = 0;
    
    if (nexFile.markers{i}.NMarkers > 0)
        % check the first marker field
        if(isfield(nexFile.markers{i}.values{1,1}, 'numericValues'))
            % convert to strings
            disp('markers are stored as numbers, converting marker values to strings');
            for field = 1:nexFile.markers{i}.NMarkers
                for item = 1:size(nexFile.markers{i}.values{field,1}.numericValues, 1)
                    nexFile.markers{i}.values{field,1}.strings{item,1} = sprintf('%d', nexFile.markers{i}.values{field,1}.numericValues(item));
                end
            end
        end
    end
    
    MarkerLength = 0;
    for j = 1:nexFile.markers{i}.NMarkers
        for k = 1:size(nexFile.markers{i}.values{j,1}.strings, 1)
            MarkerLength = max(MarkerLength, size(nexFile.markers{i}.values{j,1}.strings{k,1}, 2));
        end
    end
    
    % add extra char to hold zero (end of string)
    MarkerLength = MarkerLength + 1;
    nexFile.markers{i}.MarkerLength = MarkerLength;
    
    % marker variable type is 6
    varHeader.Type = 6;
    varHeader.Version = 100;
    varHeader.Name = nexFile.markers{i}.name;
    varHeader.Count = size(nexFile.markers{i}.timestamps,1);
    varHeader.DataOffset = dataOffset;
    varHeader.NMarkers = nexFile.markers{i}.NMarkers;
    varHeader.MarkerLength = nexFile.markers{i}.MarkerLength;
    
    writeNexVarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*4 + nexFile.markers{i}.NMarkers*64 + nexFile.markers{i}.NMarkers*varHeader.Count*MarkerLength;
    if dataOffset >= 2147483647
        fclose(fid);
        delete(fileName);
        error 'Data will not fit into .nex file (2 GB limit). Save data as .nex5 file instead.'
    end
end

for i = 1:neuronCount
    fwrite(fid, nexFile.neurons{i}.timestamps.*nexFile.freq, 'int32');
end

for i = 1:eventCount
    fwrite(fid, nexFile.events{i}.timestamps.*nexFile.freq, 'int32');
end

for i = 1:intervalCount
    fwrite(fid, nexFile.intervals{i}.intStarts.*nexFile.freq, 'int32');
    fwrite(fid, nexFile.intervals{i}.intEnds.*nexFile.freq, 'int32');
end
for i = 1:waveCount
    fwrite(fid, nexFile.waves{i}.timestamps.*nexFile.freq, 'int32');
    wf = int16(nexFile.waves{i}.waveforms./nexFile.waves{i}.ADtoMV);
    fwrite(fid, wf, 'int16');
end

for i = 1:contCount
    fwrite(fid, nexFile.contvars{i}.timestamps.*nexFile.freq, 'int32');
    fwrite(fid, nexFile.contvars{i}.fragmentStarts - 1, 'int32');
    fwrite(fid, int16(nexFile.contvars{i}.data./nexFile.contvars{i}.ADtoMV), 'int16');
end

for i = 1:markerCount
    fwrite(fid, nexFile.markers{i}.timestamps.*nexFile.freq, 'int32');
    for j = 1:nexFile.markers{i}.NMarkers
        writeStringPaddedWithZeros(fid, nexFile.markers{i}.values{j,1}.name, 64);
        for k = 1:size(nexFile.markers{i}.values{j,1}.strings, 1)
            writeStringPaddedWithZeros( fid, nexFile.markers{i}.values{j,1}.strings{k,1}, nexFile.markers{i}.MarkerLength );
        end
    end
end

fclose(fid);
result = 1;
