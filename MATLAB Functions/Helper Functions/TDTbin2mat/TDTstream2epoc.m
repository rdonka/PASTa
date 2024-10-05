function data = TDTstream2epoc(data, STREAM, THRESHOLD)
%TDTSTREAM2EPOC  stream event to epoc event threshold converter
%   data = TDTstream2epoc(DATA, STREAM, THRESHOLD), where DATA is the
%   output of TDTbin2mat, STREAM is the name of the stream store to
%   convert, THRESHOLD is the test level to convert it to a logic signal.
%
%   If data.streams.(STREAM).data is multi-channel (nrows > 1), one epoc
%   for each channel is created, like 'Wav1_1', 'Wav1_2', etc.
%
%   data    contains new EPOC store(s) with onsets and offsets from logic
%           conversion.
%
%   Example
%      data = TDTbin2mat('C:\TDT\OpenEx\Tanks\DEMOTANK2\Block-1');
%      data = TDTstream2epoc(data, 'Wav1', 0.5);
%      data.epocs.Wav1.onset
%

% iterate over all chunks in stream
tr = size(data.time_ranges, 2);

% iterate over all channels in stream
if tr > 1
    nc = size(data.streams.(STREAM).data{1}, 1);
else
    nc = size(data.streams.(STREAM).data, 1);
end

FS = data.streams.(STREAM).fs;

for ttt = 1:tr
    for ii = 1:nc
        % find onset/offset samples
        if tr > 1
            ind = data.streams.(STREAM).data{ttt}(ii, :) > THRESHOLD;
        else
            ind = data.streams.(STREAM).data(ii, :) > THRESHOLD;
        end
        rise = find(diff(ind, 1, 2) > 0)';
        fall = find(diff(ind, 1, 2) < 0)';

        % set epoc name
        if nc > 1
            name = [STREAM '_' num2str(ii)];
        else
            name = STREAM;
        end

        % convert onset/offsets to time
        st = (ceil(data.time_ranges(ttt, 1) * FS) / FS);
        if ~isfield(data.epocs, name)
            data.epocs.(name).onset = rise / FS + st;
            data.epocs.(name).offset = fall / FS + st;
        else
            data.epocs.(name).onset = [data.epocs.(name).onset; rise / FS + st];
            data.epocs.(name).offset = [data.epocs.(name).offset; fall / FS + st];
        end
        % clean up
        data.epocs.(name).onset = round((data.epocs.(name).onset/2.56e-6)) * 2.56e-6;
        data.epocs.(name).offset = round((data.epocs.(name).offset/2.56e-6)) * 2.56e-6;
        if length(data.epocs.(name).onset) > length(data.epocs.(name).offset)
            data.epocs.(name).offset = [data.epocs.(name).offset; data.time_ranges(ttt, 2)];
        end
        data.epocs.(name).name = name;
        data.epocs.(name).data = ones(size(data.epocs.(name).onset));
    end
end