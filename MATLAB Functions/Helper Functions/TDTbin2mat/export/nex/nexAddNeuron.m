function [ nexFile ] = nexAddNeuron( nexFile, timestamps, name, wireNumber, unitNumber )
% [nexFile] = nexAddNeuron( nexFile, timestamps, name ) -- adds a neuron 
%             to nexFile data structure
%
% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   timestamps - vector of neuron timestamps in seconds
%   name - neuron name
%   wireNumber - wire number (optional)
%   unitNumber - unit number (optional). Zero unit number means unsorted.
    
    if isa(name, 'char') == 0
        error 'name should be of type char'
    end
    
    wire = 0;
    unit = 0;
    if nargin > 3
        wire = wireNumber;
    end
    if nargin > 4
        unit = unitNumber;
    end
    

    neuronCount = 0;
    if(isfield(nexFile, 'neurons'))
        neuronCount = size(nexFile.neurons, 1);
    end
    neuronCount = neuronCount+1;
    nexFile.neurons{neuronCount,1}.name = name;
    nexFile.neurons{neuronCount,1}.varVersion = 100;
    nexFile.neurons{neuronCount,1}.wireNumber = wire;
    nexFile.neurons{neuronCount,1}.unitNumber = unit;
    nexFile.neurons{neuronCount,1}.xPos = 0;
    nexFile.neurons{neuronCount,1}.yPos = 0;
    % timestamps should be a vector
    if size(timestamps,1) == 1
        % if row, transpose to vector
        nexFile.neurons{neuronCount,1}.timestamps = timestamps';
    else
         nexFile.neurons{neuronCount,1}.timestamps = timestamps;
    end
    % modify end of file timestamp value in file header
    if size(nexFile.neurons{neuronCount,1}.timestamps,1) > 0
        nexFile.tend = max(nexFile.tend, timestamps(end));
    end
end
