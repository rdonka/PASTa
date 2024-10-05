function [result] = writeNex5File(nexFile, fileName, writeContVarsAsFloats)
% [result] = writeNex5File(nexFile, fileName) -- write nexFile structure
% to the specified .nex5 file. returns 1 if succeeded, 0 if failed.
% 
% INPUT:
%   nexFile - a structure containing .nex file data
%
%           SOME FIELDS OF THIS STRUCTURE ARE NOT DESCRIBED
%           BELOW. IT IS RECOMMENDED THAT YOU READ A VALID .NEX FILE 
%           TO FILL THIS STRUCTURE, THEN MODIFY THE STRUCTURE AND SAVE IT.
%
%           IF YOU WANT TO CREATE NEW .NEX FILE, USE nexCreateFileData.m,
%           nexAddContinuous.m etc. See exampleSaveDataInNexFile.m.
%           
%   fileName - if empty string, will use File Save dialog
%   writeContVarsAsFloats - if 1, values for continuous and waveform 
%              variables are written
%              as floating point numbers. This option preserves data
%              precision, but increases file size.
%
%   nexFile - a structure containing .nex file data
%   nexFile.version - file version
%   nexFile.comment - file comment
%   nexFile.tbeg - beginning of recording session (in seconds)
%   
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
%               markers{i}.values{j}.strings - array of marker value
%                     strings (if values are stored as strings in the file)
%               markers{i}.values{j}.numericValues - numeric marker values
%                     (if values are stored as numbers in the file)
%

result = 0;
if nargin < 1
    error 'provide nexFile scructure as first parameter to writeNex5File'
end

if (nargin < 2 || isempty(fileName))
   [fname, pathname] = uiputfile('*.nex5', 'Save as');
    if isequal(fname,0)
     error 'File name was not selected'
   end
   fileName = fullfile(pathname, fname);
end

if nargin < 3
    writeContVarsAsFloats = 0;
end

% note 'l' option when opening the file. 
% this options means that the file is 'little-endian'.
% this should ensure that the files are written correctly 
% on big-endian systems, such as Mac G5.
fid = fopen(fileName, 'w', 'l', 'US-ASCII');
if(fid == -1)
   error 'Unable to open file'
end

% count all the variables
neuronCount = 0;
eventCount = 0;
intervalCount = 0;
waveCount = 0;
contCount = 0;
markerCount = 0;
maxTimestampTick = maxTimestampInTicks(nexFile);

% check if we can save timestamps as 32-bit integers
if maxTimestampTick > pow2(31)
    timestampDataType = 1;
    bytesInTimestamp = 8;
else
    timestampDataType = 0;
    bytesInTimestamp = 4;
end

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
fwrite(fid, 894977358, 'int32');

% set file version to 501 or 502 depending on timestamp representation
if timestampDataType == 0
    fwrite(fid, 501, 'int32'); 
else
    fwrite(fid, 502, 'int32'); 
end
    
writeStringPaddedWithZeros(fid, nexFile.comment, 256);
fwrite(fid, nexFile.freq, 'double');
fwrite(fid, int64(nexFile.tbeg*nexFile.freq), 'int64');
fwrite(fid, nvar, 'int32');
% write zero meta offset
fwrite(fid, 0, 'uint64');
% write file time end
fwrite(fid, int64(maxTimestampTick), 'int64');

% skip padding
fwrite(fid, char(zeros(1, 56)), 'char');

% calculate where variable data starts
try
    int64(5) + int64(1);
    dataOffset = int64(356 + nvar*244);
catch
    dataOffset = double(356 + nvar*244);
end

% write variable headers
varHeader.Type = 0;
varHeader.Version = 500;
varHeader.Name = ' ';
varHeader.DataOffset = 0;
varHeader.Count = int64(0);    
varHeader.TimestampDataType = timestampDataType;
varHeader.ContinuousDataType = 0;
varHeader.SamplingFrequency = 0;
varHeader.Units = 'mV';
varHeader.ADtoUnitsCoefficient = 0;
varHeader.UnitsOffset  = 0;
varHeader.NumberOfDataPoints = 0;
varHeader.PrethresholdTimeInSeconds = 0;
varHeader.MarkerDataType = 0;
varHeader.NumberOfMarkerFields = 0;
varHeader.MarkerLength = 0;
varHeader.ContinuousIndexOfFirstPointInFramgmentDataType = 0;

meta = '{"file":{"writerSoftware":{"name":"Matlab","version":"Mar-04-16"}}, "variables":[';

% write neuron headers
for i = 1:neuronCount
    varHeader.Type = 0;
    varHeader.Name = nexFile.neurons{i}.name;
    varHeader.Count = size(nexFile.neurons{i}.timestamps,1);
    varHeader.DataOffset = int64(dataOffset);
    
    writeNex5VarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*bytesInTimestamp;
    
    if(isfield(nexFile.neurons{i},'wireNumber') && isfield(nexFile.neurons{i}, 'unitNumber')) 
        s = sprintf('{"name":"%s","probe":{"wireNumber":%d},"unitNumber":%d},',nexFile.neurons{i}.name, nexFile.neurons{i}.wireNumber,nexFile.neurons{i}.unitNumber);
        meta = strcat(meta,s);
    end
end
   
% event headers
for i = 1:eventCount
    varHeader.Type = 1;
    varHeader.Name = nexFile.events{i}.name;
    varHeader.Count = size(nexFile.events{i}.timestamps,1);
    varHeader.DataOffset = int64(dataOffset);
    
    writeNex5VarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*bytesInTimestamp;
end
    
% interval headers
for i = 1:intervalCount
    varHeader.Type = 2;
    varHeader.Name = nexFile.intervals{i}.name;
    varHeader.Count = size(nexFile.intervals{i}.intStarts,1);
    varHeader.DataOffset = int64(dataOffset);
    
    writeNex5VarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*bytesInTimestamp*2;
end

% wave headers
varHeader.TimestampDataType = timestampDataType;
varHeader.ContinuousDataType = 0;
varHeader.SamplingFrequency = 0;
varHeader.Units = 'mV';
varHeader.ADtoUnitsCoefficient = 0;
varHeader.UnitsOffset  = 0;
varHeader.NumberOfDataPoints = 0;
varHeader.PrethresholdTimeInSeconds = 0;
varHeader.MarkerDataType = 0;
varHeader.NumberOfMarkerFields = 0;
varHeader.MarkerLength = 0;
varHeader.ContinuousIndexOfFirstPointInFramgmentDataType = 0;

for i = 1:waveCount
    varHeader.Type = 3;
    varHeader.Name = nexFile.waves{i}.name;
    varHeader.Count = size(nexFile.waves{i}.timestamps,1);
    varHeader.DataOffset = int64(dataOffset);
    varHeader.SamplingFrequency = nexFile.waves{i}.WFrequency;
    varHeader.NumberOfDataPoints = nexFile.waves{i}.NPointsWave;
    varHeader.PrethresholdTimeInSeconds = 0;
    if(isfield(nexFile.waves{i}, 'PrethresholdTimeInSeconds'))
        varHeader.PrethresholdTimeInSeconds = nexFile.waves{i}.PrethresholdTimeInSeconds;
    end
    
    dataPointSize = 2;
    
    if writeContVarsAsFloats == 1
        varHeader.ContinuousDataType = 1;
        nexFile.waves{i}.ADtoMV = 1;
        dataPointSize = 4;
    else
        varHeader.ContinuousDataType = 0;
        c = 1;
        if varHeader.Count > 0
            % we need to recalculate a/d to milliVolts factor
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
    end
    
    nexFile.waves{i}.MVOfffset = 0;
    
    varHeader.ADtoUnitsCoefficient = nexFile.waves{i}.ADtoMV;
    varHeader.UnitsOffset = nexFile.waves{i}.MVOfffset;
    
    writeNex5VarHeader(fid, varHeader);
    
    dataOffset = dataOffset + varHeader.Count*bytesInTimestamp + varHeader.NumberOfDataPoints*varHeader.Count*dataPointSize;
    
    if(isfield(nexFile.waves{i},'wireNumber') && isfield(nexFile.waves{i}, 'unitNumber')) 
        s = sprintf('{"name":"%s","probe":{"wireNumber":%d},"unitNumber":%d},',nexFile.waves{i}.name, nexFile.waves{i}.wireNumber,nexFile.waves{i}.unitNumber);
        meta = strcat(meta,s);
    end
end

% continuous variables
varHeader.TimestampDataType = timestampDataType;
varHeader.ContinuousDataType = 0;
varHeader.SamplingFrequency = 0;
varHeader.Units = 'mV';
varHeader.ADtoUnitsCoefficient = 0;
varHeader.UnitsOffset  = 0;
varHeader.NumberOfDataPoints = 0;
varHeader.PrethresholdTimeInSeconds = 0;
varHeader.MarkerDataType = 0;
varHeader.NumberOfMarkerFields = 0;
varHeader.MarkerLength = 0;
varHeader.ContinuousIndexOfFirstPointInFramgmentDataType = 0;

for i = 1:contCount
    varHeader.Type = 5;
    varHeader.Name = nexFile.contvars{i}.name;
    varHeader.Count = size(nexFile.contvars{i}.timestamps,1);
    varHeader.DataOffset = int64(dataOffset);
    varHeader.SamplingFrequency = nexFile.contvars{i}.ADFrequency;
    varHeader.NumberOfDataPoints = size(nexFile.contvars{i}.data, 1);
    nexFile.contvars{i}.MVOfffset = 0;
    varHeader.UnitsOffset = nexFile.contvars{i}.MVOfffset;
    dataPointSize = 2;
    
    if writeContVarsAsFloats == 1
        varHeader.ContinuousDataType = 1;
        nexFile.contvars{i}.ADtoMV = 1;
        dataPointSize = 4;
    else
        varHeader.ContinuousDataType = 0;
        c = 1;
        if varHeader.Count > 0
            % we need to recalculate a/d to milliVolts factor
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
    end
    
    varHeader.ADtoUnitsCoefficient = nexFile.contvars{i}.ADtoMV;        

    writeNex5VarHeader(fid, varHeader);

    dataOffset = dataOffset + varHeader.Count*(bytesInTimestamp+4) + varHeader.NumberOfDataPoints*dataPointSize;
end

% markers
varHeader.TimestampDataType = timestampDataType;
varHeader.ContinuousDataType = 0;
varHeader.SamplingFrequency = 0;
varHeader.Units = 'mV';
varHeader.ADtoUnitsCoefficient = 0;
varHeader.UnitsOffset  = 0;
varHeader.NumberOfDataPoints = 0;
varHeader.PrethresholdTimeInSeconds = 0;
varHeader.MarkerDataType = 0;
varHeader.NumberOfMarkerFields = 0;
varHeader.MarkerLength = 0;
varHeader.ContinuousIndexOfFirstPointInFramgmentDataType = 0;

for i = 1:markerCount
    varHeader.MarkerDataType = 0;
    
    nexFile.markers{i}.NMarkers = size(nexFile.markers{i}.values, 1);
    nexFile.markers{i}.MarkerLength = 0;  
    if (nexFile.markers{i}.NMarkers > 0)
        % check the first marker field
        if(isfield(nexFile.markers{i}.values{1,1}, 'numericValues')) 
            varHeader.MarkerDataType = 1;
        end
    end
    
    if (varHeader.MarkerDataType == 0)
        MarkerLength = 0;
        for j = 1:nexFile.markers{i}.NMarkers
          for k = 1:size(nexFile.markers{i}.values{j,1}.strings, 1)
            MarkerLength = max(MarkerLength, size(nexFile.markers{i}.values{j,1}.strings{k,1}, 2));
          end
        end
        % add extra char to hold zero (end of string)
        MarkerLength = MarkerLength + 1;
        nexFile.markers{i}.MarkerLength = MarkerLength;    
    end    
    
    varHeader.Type = 6;
    varHeader.Name = nexFile.markers{i}.name;
    varHeader.Count = size(nexFile.markers{i}.timestamps,1);
    varHeader.DataOffset = int64(dataOffset);
    varHeader.SamplingFrequency = 0;
    varHeader.ADtoUnitsCoefficient = 0;
    varHeader.UnitsOffset = 0;
    varHeader.NumberOfDataPoints = 0;
    varHeader.NumberOfMarkerFields = nexFile.markers{i}.NMarkers;
    varHeader.MarkerLength = nexFile.markers{i}.MarkerLength;
    
    writeNex5VarHeader(fid, varHeader);
    if (varHeader.MarkerDataType == 0)
        dataOffset = dataOffset + varHeader.Count*bytesInTimestamp + varHeader.NumberOfMarkerFields*64 + varHeader.NumberOfMarkerFields*varHeader.Count*varHeader.MarkerLength;    
    else
        % we have 4 bytes per marker value
        dataOffset = dataOffset + varHeader.Count*bytesInTimestamp + varHeader.NumberOfMarkerFields*64 + varHeader.NumberOfMarkerFields*varHeader.Count*4; 
    end
end

% write variable data after headers
if timestampDataType == 0
    timestampFormat = 'int32';
else
    timestampFormat = 'int64';
end

for i = 1:neuronCount
    fwrite(fid, nexFile.neurons{i}.timestamps.*nexFile.freq, timestampFormat);
end

for i = 1:eventCount
    fwrite(fid, nexFile.events{i}.timestamps.*nexFile.freq, timestampFormat);
end

for i = 1:intervalCount
    fwrite(fid, nexFile.intervals{i}.intStarts.*nexFile.freq, timestampFormat);
    fwrite(fid, nexFile.intervals{i}.intEnds.*nexFile.freq, timestampFormat);
end

for i = 1:waveCount
    fwrite(fid, nexFile.waves{i}.timestamps.*nexFile.freq, timestampFormat);
    if writeContVarsAsFloats == 1
        fwrite(fid, single(nexFile.waves{i}.waveforms), 'float');
    else
        fwrite(fid, int16(nexFile.waves{i}.waveforms./nexFile.waves{i}.ADtoMV), 'int16');
    end
end

for i = 1:contCount
    fwrite(fid, nexFile.contvars{i}.timestamps.*nexFile.freq, timestampFormat);
    fwrite(fid, nexFile.contvars{i}.fragmentStarts - 1, 'int32');
    if writeContVarsAsFloats == 1
        fwrite(fid, single(nexFile.contvars{i}.data), 'float');
    else
        fwrite(fid, int16(nexFile.contvars{i}.data./nexFile.contvars{i}.ADtoMV), 'int16');
    end
end

for i = 1:markerCount
    fwrite(fid, nexFile.markers{i}.timestamps.*nexFile.freq, timestampFormat);
    for j = 1:nexFile.markers{i}.NMarkers
      writeStringPaddedWithZeros(fid, nexFile.markers{i}.values{j,1}.name, 64);
      if( isfield(nexFile.markers{i}.values{j,1}, 'numericValues') ) 
          fwrite(fid, nexFile.markers{i}.values{j,1}.numericValues, 'uint32');           
      else
          for k = 1:size(nexFile.markers{i}.values{j,1}.strings, 1)
            writeStringPaddedWithZeros( fid, nexFile.markers{i}.values{j,1}.strings{k,1}, nexFile.markers{i}.MarkerLength );
          end    
      end      
    end
end

% write metadata
if meta(end) == ','
    meta(end) = ' ';
end
meta = strcat(meta,']}');
position = ftell(fid);
fwrite(fid, meta, 'char');
metaPos = 284;
fseek(fid, metaPos, 'bof');
fwrite(fid, position, 'uint64');
fclose(fid);
result = 1;
