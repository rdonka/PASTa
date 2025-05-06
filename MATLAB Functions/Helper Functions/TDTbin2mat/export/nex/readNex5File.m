function [nexFile] = readNex5File(fileName)
% [nexFile] = readNexFile(fileName) -- read .nex5 file and return file data
%             in nexFile structure
%
% INPUT:
%   fileName - if empty string, will use File Open dialog
%
% OUTPUT:
%   nexFile - a structure containing .nex file data
%   nexFile.version - file version
%   nexFile.comment - file comment
%   nexFile.freq - file timestamp frequency (Hz)
%   nexFile.tbeg - beginning of recording session (in seconds)
%   nexFile.tend - end of recording session (in seconds)
%   nexFile.metadata - file metadata as a string in json format
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
%                                waveform is a column 
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
%               markers{i}.values{j}.name - name of marker value 
%               markers{i}.values{j}.strings - array of marker value strings 
%                     (if values are stored as strings in the file)
%               markers{i}.values{j}.numericValues - numeric marker values
%                     (if values are stored as numbers in the file)
%

nexFile = [];

if (nargin == 0 | isempty(fileName))
   [fname, pathname] = uigetfile('*.nex5', 'Select a .nex5 NeuroExplorer file');
   if isequal(fname,0)
     error 'No file was selected'
     return
   end
   fileName = fullfile(pathname, fname);
end

% note 'l' option when opening the file. 
% this options means that the file is 'little-endian'.
% this should ensure that the files are read correctly 
% on big-endian systems, such as Mac G5.
fid = fopen(fileName, 'r', 'l','US-ASCII');
if(fid == -1)
   error 'Unable to open file'
   return
end

magic = fread(fid, 1, 'int32');
if magic ~= 894977358
    error 'The file is not a valid .nex5 file'
end
nexFile.nex5Version = fread(fid, 1, 'int32'); 
comment = fread(fid, 256, '*char')'; 
% remove first zero and all characters after the first zero
comment(end+1) = 0; 
nexFile.comment = comment(1:min(find(comment==0))-1);
nexFile.freq = fread(fid, 1, 'double');
nexFile.tbeg = fread(fid, 1, 'int64')./nexFile.freq;
nvar = fread(fid, 1, 'int32');
metaOffset = fread(fid, 1, 'uint64');
nexFile.tend = fread(fid, 1, 'int64')./nexFile.freq;
if nexFile.nex5Version < 501
    nexFile.tend = 0;        
end
% skip padding
fseek(fid, 56, 'cof');    

neuronCount = 0;
eventCount = 0;
intervalCount = 0;
waveCount = 0;
popCount = 0;
contCount = 0;
markerCount = 0;

% real all variables
for variableIndex=1:nvar
    % read variable header
    type = fread(fid, 1, 'int32');
    varVersion = fread(fid, 1, 'int32');
    name = fread(fid, 64, '*char')';
    % remove first zero and all characters after the first zero
    name(end+1) = 0;
    name = name(1:min(find(name==0))-1);
    offset = fread(fid, 1, 'uint64');
    n = fread(fid, 1, 'int64');
    tsDataType = fread(fid, 1, 'int32');
    contDataType = fread(fid, 1, 'int32');
    WFrequency = fread(fid, 1, 'double'); % wf sampling fr.
    units = fread(fid, 32, '*char')';
    units(end+1) = 0;
    units = units(1:min(find(units==0))-1);
    ADtoMV  = fread(fid, 1, 'double'); % coeff to convert from AD values to Millivolts.
    MVOfffset = fread(fid, 1, 'double'); % coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset
    NPointsWave = fread(fid, 1, 'uint64'); % number of points in each wave
    PreThresholdTime = fread(fid, 1, 'double'); 

    markerDataType = fread(fid, 1, 'int32'); 

    NMarkers = fread(fid, 1, 'int32'); % how many values are associated with each marker
    MarkerLength = fread(fid, 1, 'int32'); % how many characters are in each marker value
    contFragmentStartDataType = fread(fid, 1, 'int32');
    filePosition = ftell(fid);
    switch type
        case 0 % neuron
            neuronCount = neuronCount+1;
            nexFile.neurons{neuronCount,1}.name = name;
            nexFile.neurons{neuronCount,1}.varNex5Version = varVersion;
            % go to the variable data position and read timestamps
            fseek(fid, offset, 'bof');
            if ( tsDataType == 0 )
                nexFile.neurons{neuronCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            else
                nexFile.neurons{neuronCount,1}.timestamps = fread(fid, [n 1], 'int64')./nexFile.freq;
            end
                        
        case 1 % event
            eventCount = eventCount+1;
            nexFile.events{eventCount,1}.name = name;
            nexFile.events{eventCount,1}.varNex5Version = varVersion;
            fseek(fid, offset, 'bof');
            if ( tsDataType == 0 )
                nexFile.events{eventCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            else
                nexFile.events{eventCount,1}.timestamps = fread(fid, [n 1], 'int64')./nexFile.freq;
            end
                    
        case 2 % interval
            intervalCount = intervalCount+1;
            nexFile.intervals{intervalCount,1}.name = name;
            nexFile.intervals{intervalCount,1}.varNex5Version = varVersion;
            fseek(fid, offset, 'bof');
            if ( tsDataType == 0 )
                nexFile.intervals{intervalCount,1}.intStarts = fread(fid, [n 1], 'int32')./nexFile.freq;
                nexFile.intervals{intervalCount,1}.intEnds = fread(fid, [n 1], 'int32')./nexFile.freq;
            else
                nexFile.intervals{intervalCount,1}.intStarts = fread(fid, [n 1], 'int64')./nexFile.freq;
                nexFile.intervals{intervalCount,1}.intEnds = fread(fid, [n 1], 'int64')./nexFile.freq;
            end
                    
        case 3 % waveform
            waveCount = waveCount+1;
            nexFile.waves{waveCount,1}.name = name;
            nexFile.waves{waveCount,1}.varNex5Version = varVersion;
            nexFile.waves{waveCount,1}.NPointsWave = NPointsWave;
            nexFile.waves{waveCount,1}.WFrequency = WFrequency;
            nexFile.waves{waveCount,1}.ADtoMV = ADtoMV;
            nexFile.waves{waveCount,1}.MVOfffset = MVOfffset;
            nexFile.waves{waveCount,1}.PreThresholdTime = PreThresholdTime;
            fseek(fid, offset, 'bof');
            if ( tsDataType == 0 )
                nexFile.waves{waveCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            else
                nexFile.waves{waveCount,1}.timestamps = fread(fid, [n 1], 'int64')./nexFile.freq;
            end
            if ( contDataType == 0 )
                wf = fread(fid, [NPointsWave n], 'int16');
            else
                wf = fread(fid, [NPointsWave n], 'float32');
            end
            nexFile.waves{waveCount,1}.waveforms = wf.*ADtoMV + nexFile.waves{waveCount,1}.MVOfffset;
                        
        case 4 % population vector
            popCount = popCount+1;
            nexFile.popvectors{popCount,1}.name = name;
            nexFile.popvectors{popCount,1}.varNex5Version = varVersion;
            fseek(fid, offset, 'bof');
            nexFile.popvectors{popCount,1}.weights = fread(fid, [n 1], 'double');
                        
        case 5 % continuous variable
            contCount = contCount+1;
            nexFile.contvars{contCount,1}.name = name;
            nexFile.contvars{contCount,1}.varNex5Version = varVersion;
            nexFile.contvars{contCount,1}.ADtoMV = ADtoMV;
            nexFile.contvars{contCount,1}.MVOfffset = MVOfffset;
            nexFile.contvars{contCount,1}.ADFrequency = WFrequency;
            fseek(fid, offset, 'bof');
            if ( tsDataType == 0 )
                nexFile.contvars{contCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            else
                nexFile.contvars{contCount,1}.timestamps = fread(fid, [n 1], 'int64')./nexFile.freq;
            end
            if ( contFragmentStartDataType == 0 )
                nexFile.contvars{contCount,1}.fragmentStarts = fread(fid, [n 1], 'uint32') + 1;
            else
                nexFile.contvars{contCount,1}.fragmentStarts = fread(fid, [n 1], 'uint64') + 1;
            end
            
            if ( contDataType == 0 )
                nexFile.contvars{contCount,1}.data = fread(fid, [NPointsWave 1], 'int16').*ADtoMV + nexFile.contvars{contCount,1}.MVOfffset;
            else
                nexFile.contvars{contCount,1}.data = fread(fid, [NPointsWave 1], 'float32').*ADtoMV + nexFile.contvars{contCount,1}.MVOfffset;
            end
                     
        case 6 % marker
            markerCount = markerCount+1;
            nexFile.markers{markerCount,1}.name = name;
            nexFile.markers{markerCount,1}.varNex5Version = varVersion;
            fseek(fid, offset, 'bof');
            if ( tsDataType == 0 )
                nexFile.markers{markerCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
            else
                nexFile.markers{markerCount,1}.timestamps = fread(fid, [n 1], 'int64')./nexFile.freq;
            end
            for markerFieldIndex=1:NMarkers
                markerName = fread(fid, 64, '*char')';
                markerName(end+1) = 0;
                markerName = markerName(1:min(find(markerName==0))-1);
                nexFile.markers{markerCount,1}.values{markerFieldIndex,1}.name = markerName;
                if ( markerDataType == 0 )
                    for markerValueIndex = 1:n
                        markerValue = fread(fid, MarkerLength, '*char')';
                        % remove first zero and all characters after the first zero
                        markerValue(end+1) = 0;
                        markerValue = markerValue(1:min(find(markerValue==0))-1);
                        nexFile.markers{markerCount,1}.values{markerFieldIndex,1}.strings{markerValueIndex, 1} = markerValue;
                    end
                else
                    nexFile.markers{markerCount,1}.values{markerFieldIndex,1}.numericValues = fread(fid, [n 1], 'uint32');
                end
            end
                    
        otherwise
            disp (['unknown variable type ' num2str(type)]);
    end
    % return to file position that was after reading the variable header
    fseek(fid, filePosition, 'bof');
    dummy = fread(fid, 60, 'char');
end

% read and process metadata at the end of the file
try
    fseek(fid, 0,'eof');
    filelength = ftell(fid);
    if ( metaOffset > 0 && metaOffset < filelength )
        fseek(fid, metaOffset, 'bof');
        jsonMeta = fread(fid, filelength - metaOffset, '*char')';
        jsonMeta(end+1) = 0;
        jsonMeta = jsonMeta(1:min(find(jsonMeta==0))-1);  
        nexFile.metadata = jsonMeta;
        % disp(jsonMeta);
        meta = parse_json(nexFile.metadata);
        meta = meta{1};
        for i=1:length(meta.variables)
            varMeta = meta.variables{i};
            name = varMeta.name;
            if isfield(nexFile, 'neurons')
                for j=1:length(nexFile.neurons)
                    nrName = nexFile.neurons{j}.name;
                    if strcmp(nrName,name) == 1
                        if isfield(varMeta, 'unitNumber')
                            nexFile.neurons{j}.unitNumber = varMeta.unitNumber;
                        end
                        
                        if isfield(varMeta, 'probe')
                            if isfield(varMeta.probe, 'wireNumber')
                                nexFile.neurons{j}.wireNumber = varMeta.probe.wireNumber;
                            end
                            if isfield(varMeta.probe, 'position')
                                nexFile.neurons{j}.xPos = varMeta.probe.position.x;
                                nexFile.neurons{j}.yPos = varMeta.probe.position.y;
                            end
                        end                       
                        break
                    end
                end 
            end
            if isfield(nexFile, 'waves')
                for j=1:length(nexFile.waves)
                    waveName = nexFile.waves{j}.name;
                    if strcmp(waveName,name) == 1
                        if isfield(varMeta, 'unitNumber')
                            nexFile.waves{j}.unitNumber = varMeta.unitNumber;
                        end
                        if isfield(varMeta, 'probe')
                            if isfield(varMeta.probe, 'wireNumber')
                                nexFile.waves{j}.wireNumber = varMeta.probe.wireNumber;
                            end
                        end
                        break
                    end
                end  
            end
        end
    end
catch ME
    msgText = getReport(ME)
    warning('unable to process metadata')
    warning(msgText)
end
fclose(fid);