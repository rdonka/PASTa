function [ maxTimestamp ] = maxTimestampInTicks( nexFile )
%[result] = maxTimestampInTicks(nedFile) -- calculates maximum timestamp
% of the variables in nexFile data structure.

maxTimestamp = 0;

if(isfield(nexFile, 'neurons'))
    for i = 1:size(nexFile.neurons, 1)
        if size(nexFile.neurons{i}.timestamps,1) > 0
            maxTimestamp = max(maxTimestamp, nexFile.neurons{i}.timestamps(end));
        end
    end
end

if(isfield(nexFile, 'events'))
    for i = 1:size(nexFile.events, 1)
        if size(nexFile.events{i}.timestamps,1) > 0
            maxTimestamp = max(maxTimestamp, nexFile.events{i}.timestamps(end));
        end
    end
end

if(isfield(nexFile, 'intervals'))
    for i = 1:size(nexFile.intervals, 1)
        if size(nexFile.intervals{i}.intEnds,1) > 0
            maxTimestamp = max(maxTimestamp, nexFile.intervals{i}.intEnds(end));
        end
    end
end

if(isfield(nexFile, 'waves'))
    for i = 1:size(nexFile.waves, 1)
        if size(nexFile.waves{i}.timestamps,1) > 0
            step = 1.0/nexFile.waves{i}.WFrequency;
            maxTimestamp = max(maxTimestamp, nexFile.waves{i}.timestamps(end)+step*(nexFile.waves{i}.NPointsWave-1));
        end
    end
end

if(isfield(nexFile, 'contvars'))
    for i = 1:size(nexFile.contvars, 1)
        if size(nexFile.contvars{i}.timestamps,1) > 0
            step = 1.0/nexFile.contvars{i}.ADFrequency;
            numPointsInLastFragment = size(nexFile.contvars{i}.data, 1) - nexFile.contvars{i}.fragmentStarts(end) - 1;
            maxTimestamp = max(maxTimestamp, nexFile.contvars{i}.timestamps(end)+step*(numPointsInLastFragment-1));
        end
    end
end

if(isfield(nexFile, 'markers'))
    for i = 1:size(nexFile.markers, 1)
        if size(nexFile.markers{i}.timestamps,1) > 0
            maxTimestamp = max(maxTimestamp, nexFile.markers{i}.timestamps(end));
        end
    end
end

maxTimestamp = maxTimestamp  * nexFile.freq;

end

