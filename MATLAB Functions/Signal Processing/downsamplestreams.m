function [data] = downsamplestreams(data,streamfieldnames,fsfieldname,targetfs,varargin)
% FINDSUBJECTMEANTRIALTRACE    Finds the mean trace by trial type for each subject.
%
% INPUTS:
%   DATA:               This is a structure that contains at least the
%                           specified field with trial by trial signal and 
%                           trial type.
%
%   STREAMFIELDNAMES:   Cell array containing the names of the field in the data structure that
%                       contains the streams to be downsampled.
%       
%   FSFIELDNAME:  The name of the field in the data structure
%                            that contains the sampling rate of the data
%                            streams.
%
%   TARGETFS: Cell array; names of fields containing
%                               relevant trial epochs. These will be averaged and rounded to
%                           identify common mean event epochs.
%
% OPTIONAL INPUTS:
%   DSMETHOD: String; Output name for the trial ID field
%                          output with the mean trace data. Default:
%                          "TrialType".
%
% OUTPUTS:
%       DATADS:           This is the original data structure with added 
%                       fields for the mean signal for each trial type.
%                       Fields will be named as 
%                       data.mean_TRIALSTREAMFIELDNAME with each trial type in the field "TrialType".
% 
% Written by R M Donka, January 2024
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

    %% Prepare Settings
    % Import required and optional inputs into a structure
    defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'dsmethod', defaultparameters.dsmethod); 
    addParameter(p, 'adjustepocfieldnames', []); 

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Initialize params
    params.streamfieldnames = streamfieldnames;
    params.fsfieldname = fsfieldname;
    params.targetfs = targetfs;

    dsmethod = params.dsmethod;
    adjustepocfieldnames = params.adjustepocfieldnames;

    % Prepare output ds field name    
    disp(['DOWNSAMPLESTREAM: Downsampling raw data streams to target fs of',num2str(targetfs),'.'])
    disp(params)

    for eachfile = 1:length(data)
        disp(['Downsampling stream: File ',num2str(eachfile)]) % Display which file is being processed

        startingfs = data(eachfile).(fsfieldname);
        dsFactor = round(startingfs / targetfs); % The decimation factor
        
        for eachstream = 1:length(streamfieldnames)
            streamfieldname = streamfieldnames{eachstream};
            currrawstreamdata = data(eachfile).(streamfieldname);
    
            startingnsamples = length(currrawstreamdata);
            dsnsamples = floor(startingnsamples/dsFactor);
    
            trimrawstreamdata = currrawstreamdata(:, 1:(dsnsamples*dsFactor));
            reshapedrawstreamdata = reshape(trimrawstreamdata, dsFactor, dsnsamples);
    
            if strcmp(dsmethod, 'mean')
                dsstreamdata = mean(reshapedrawstreamdata,1);
            elseif strcmp(dsmethod, 'median')
                dsstreamdata = median(reshapedrawstreamdata,1);
            else
                error('ERROR: dsmethod not recognized. Options are "mean" and "median".')
            end
    
            data(eachfile).(streamfieldname) = dsstreamdata;
        end

        if ~isempty(adjustepocfieldnames)
            for eachepoc = 1:length(adjustepocfieldnames)
                currepocfieldname = adjustepocfieldnames{eachepoc};
    
                if isfield(data(eachfile), currepocfieldname)
                    rawidxs = data(eachfile).(currepocfieldname);
                    
                    if ~isempty(rawidxs)
                        % Scale indices
                        data(eachfile).(currepocfieldname) = floor(rawidxs ./ dsFactor);
                    end
                end
            end
        end

        data(eachfile).fs = targetfs;
        data(eachfile).params.fs = targetfs;
    end
end