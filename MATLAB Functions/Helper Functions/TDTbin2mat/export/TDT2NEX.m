close all; clear all; clc; % cleanup the workspace
addpath('nex') % add child directory for NEX export functions in it
addpath('../') % add parent directory with TDTbin2mat in it



% change to your block path
BLOCKPATH = 'C:\\TDT\\TDTExampleData\\Algernon-180308-130351';

% set to 0 to only export timestamps (as 'Neurons') in NEX
EXPORT_SNIPPET_WAVEFORMS = 1; 

% set to 1 to organize output by sort code
EXPORT_SNIPPET_SORTCODES = 0;

% read TDT data into Matlab - add any data filters here
data = TDTbin2mat(BLOCKPATH);




% BEGIN EXPORT SCRIPT
[pathstr,name,ext] = fileparts(BLOCKPATH);
nex5FilePath = strcat(BLOCKPATH, '\', name, '.nex5');
fprintf('exporting %s to %s\n', BLOCKPATH, nex5FilePath);

% start new nex file data
nexFile = nexCreateFileData(24414.0625);

% add streams
stores = fieldnames(data.streams);
for store = 1:length(stores)
    this_store = data.streams.(stores{store});
    fprintf('adding %s\n', stores{store});
    for ch = 1:size(this_store.data,1)
        varname = sprintf('%s_%03d', this_store.name, ch);
        fprintf('\tchannel %d\n', ch);
        if isa(this_store.data(ch, 1), 'int16')
            scale = 1e-3; % assume scale factor was 1e6, convert to mV
        else
            scale = 1e3; % convert V to mV
        end
        nexFile = nexAddContinuous(nexFile, ...
            data.time_ranges(1) + 1/this_store.fs, ...
            this_store.fs, ...
            single(this_store.data(ch, :)) * scale, ...
            varname);
    end
end

% add epocs
if ~isempty(data.epocs)
    stores = fieldnames(data.epocs);
    for store = 1:length(stores)
        this_store = data.epocs.(stores{store});
        fprintf('adding %s\n', this_store.name);
        % convert data values to strings
        str_data = arrayfun(@num2str, this_store.data, 'UniformOutput', 0);
        nexFile = nexAddMarker(nexFile, ...
            this_store.onset, ...
            this_store.name, ...
            str_data, ...
            this_store.name);
        nexFile = nexAddInterval(nexFile, ...
            this_store.onset, ...
            this_store.offset, ...
            ['int_' this_store.name]);
    end
end

% add snippets
if ~isempty(data.snips)
    stores = fieldnames(data.snips);
    for store = 1:length(stores)
        this_store = data.snips.(stores{store});
        fprintf('adding %s\n', this_store.name);
        numberOfPointsInWaveform = length(data.snips.(stores{store}).data(1,:));
        preThresholdTimeInSeconds = round(1/4 * numberOfPointsInWaveform) / data.snips.(stores{store}).fs;
        unique_ch = sort(unique(data.snips.(stores{store}).chan));
        for ch = 1:length(unique_ch)
            ind = data.snips.eNe1.chan == unique_ch(ch);
            fprintf('\tchannel %d\n', unique_ch(ch));
            if EXPORT_SNIPPET_SORTCODES
                unique_sc = sort(unique(data.snips.(stores{store}).sortcode(ind)));
                for sc = 1:length(unique_sc)
                    ind2 = ind & (data.snips.(stores{store}).sortcode == unique_sc(sc));
                    ts = data.snips.(stores{store}).ts(ind2);
                
                    fprintf('\t\tsortcode %d\n', unique_sc(sc));
                    name = sprintf('%s_%03d_%02d',stores{store},unique_ch(ch), unique_sc(sc));
                    if EXPORT_SNIPPET_WAVEFORMS
                        nexFile = nexAddWaveform(nexFile, ...
                            data.snips.(stores{store}).fs, ...
                            ts, ...
                            1e3 * data.snips.(stores{store}).data(ind2,:)', ...
                            name, ...
                            preThresholdTimeInSeconds, ...
                            numberOfPointsInWaveform, ...
                            unique_ch(ch), ...
                            unique_sc(sc));
                    else
                        nexFile = nexAddNeuron(nexFile, ...
                            ts, ...
                            name, ...
                            unique_ch(ch), ...
                            unique_sc(sc));
                    end
                end
            else
                ts = data.snips.(stores{store}).ts(ind);
                name = sprintf('%s_%03d',stores{store},unique_ch(ch));
                if EXPORT_SNIPPET_WAVEFORMS
                    nexFile = nexAddWaveform(nexFile, ...
                        data.snips.(stores{store}).fs, ...
                        ts, ...
                        1e3 * data.snips.(stores{store}).data(ind,:)', ...
                        name, ...
                        preThresholdTimeInSeconds, ...
                        numberOfPointsInWaveform, ...
                        unique_ch(ch));
                else
                    nexFile = nexAddNeuron(nexFile, ...
                        ts, ...
                        name, ...
                        unique_ch(ch));
                end
            end
        end
    end
end

writeNex5File(nexFile, nex5FilePath);
fprintf('writing %s\n', nex5FilePath);
disp('done')