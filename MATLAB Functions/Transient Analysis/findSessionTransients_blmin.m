function [data] = findSessionTransients_blmin(data,whichstream,whichthreshold,whichfs,thresholdlabel,varargin)
% FINDSESSIONTRANSIENTS_BLMIN   Finds transients for the whole session. 
%                               Pre-transient baselines are set to the min
%                               value within the pre-transient window.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the data
%                       stream you want to analyze.
%
%       WHICHSTREAM:    A variable containing a string with the name of the 
%                       field containing the stream to be analyzed for 
%                       transients.
%
%       WHICHFS:        The name of the field containing the sampling rate
%                       of the streams (fs).
%
%       THRESHOLDLABEL: A variable containing a string with the name of the 
%                       field containing the threshold criterion for the 
%                       stream to be analyzed for transients.
%
% OPTIONAL INPUTS:
%       PREMINSTARTMS:  Number of millseconds pre-transient to use as the
%                       start of the baseline window.
%                       Default: 8000
%
%       PREMINENDMS:    Number of millseconds pre-transient to use as the
%                       end of the baseline window.
%                       Default: 100
%
%       POSTTRANSIENTMS: Number of millseconds post-transient to use for
%                       the post peak baseline and trimmed data output.
%                       Default: 1600
%
%       QUANTIFICATIONHEIGHT: The height at which to characterize rise time,
%                       fall time, and AUC. Must be a number between 0 and 1.
%                       Default: 0.5
%  
%       OUTPUTTRANSIENTDATA: Set to 1 to output cut data streams for each
%                       transient event. Set to 0 to skip.
%                       Default: 1
%
% OUTPUTS:
%       DATA:           The original data structure with
%                       sessiontransients_blmin_THRESHOLDLABEL added in.
%                       For more details on individual variables, see the
%                       PASTa user guide. 
%                       The output contains four nested tables: 
%                       INPUTS: Includes all required and optional inputs.
%                           If optional inputs are not specified, defaults
%                           will be applied.
%                       TRANSIENTQUANTIFICATION: Includes the quantified
%                           variables for each transient, including
%                           amplitude, rise time, fall time, width, and
%                           AUC. 
%                       TRANSIENTSTREAMLOCS: Pre-transient baseline, 
%                           transient peak, rise, and fall locations for 
%                           each transient to match the cut transient 
%                           stream data.
%                       TRANSIENTSTREAMDATA: Cut data stream from baseline
%                           start to the end of the post-transient period
%                           for each transient event.
%                       Note that for all data outputs, each transient is
%                       in a separate row. If OUTPUTTRANSIENTDATA is set to 
%                       anything other than 1, the TRANSIENTSTREAMLOCS and
%                       TRANSIENTSTREAMDATA tables will be skipped and not
%                       included in the output.
%
% Written by R M Donka, October 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

disp('SESSION TRANSIENTS: Peak baseline determined by min value in the specified baseline window.')

%% Prepare Settings
% Import required and optional inputs into a structure
    inputs = struct(...
        'whichstream',[],...
        'whichthreshold',[],...
        'whichfs',[],...
        'thresholdlabel',[],...
        'preminstartms', [],...
        'preminendms',[],...
        'posttransientms',[],...
        'quantificationheight',[],...
        'outputtransientdata',[]);
    inputs = parseArgsLite(varargin,inputs);
    
    % Prepare defaults and check for optional inputs
    inputs.whichstream = whichstream;
    inputs.whichthreshold = whichthreshold;
    inputs.whichfs = whichfs;
    inputs.thresholdlabel = thresholdlabel;

    if isempty(inputs.preminstartms) % Pre-transient baseline ms start
        preminstartms = 1000;
        inputs.preminstartms = preminstartms;
    else
        preminstartms = inputs.preminstartms;
    end
    if isempty(inputs.preminendms) % Pre-transient baseline ms end
        preminendms = 100;
        inputs.preminendms = preminendms;
    else
        preminendms = inputs.preminendms;
    end
    if isempty(inputs.posttransientms) % Post transient ms for fall/AUC and output data
        posttransientms = 2000;
        inputs.posttransientms = posttransientms;
    else
        posttransientms = inputs.posttransientms;
    end
    if isempty(inputs.quantificationheight) % Height for rise/fall/AUC; defaults to 0.5 (half height)
        quantificationheight = 0.5; 
        inputs.quantificationheight = quantificationheight;
    else
        quantificationheight = inputs.quantificationheight;
    end
    if isempty(inputs.outputtransientdata)
        outputtransientdata = 1; % Height for rise/fall/AUC defaults to 0.5 (half height)
        inputs.outputtransientdata = outputtransientdata;
    else
        outputtransientdata = inputs.outputtransientdata;
    end

    disp('INPUTS:') % Display all input values
    disp(inputs)
    
    %% Find transients
    for eachfile = 1:length(data)
        disp(['Finding Transients: File ',num2str(eachfile)]) % Display which file is being processed
        
        % Prep variables
        fs = data(eachfile).(whichfs);
        blstartsamples = floor(fs*(preminstartms/1000));
        blendsamples = floor(fs*(preminendms/1000));
        posttransientsamples = floor(fs*(posttransientms/1000));

        % Find all maxes and mins in stream
        allmaxlocs = find(islocalmax(data(eachfile).(whichstream))); % Find all maxes

        % Prepare tables - preallocate size
        allvarnames = {'transientID','maxloc', 'maxval', 'preminloc', 'preminval', 'amp','risestartloc','risestartval','risesamples', 'risems','fallendloc','fallendval','fallsamples','fallms','widthsamples','widthms','AUC'};
        [allvartypes{1:length(allvarnames)}] = deal('double');
        transientquantification = table('Size',[length(allmaxlocs), length(allvarnames)], 'VariableNames', allvarnames, 'VariableTypes', allvartypes);

        if outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, prep data table and structure for individual transient streams
            transientstreamlocsvarnames = {'transientID','maxloc','preminloc','risestartloc','fallendloc'};
            [transientstreamlocsvartypes{1:length(transientstreamlocsvarnames)}] = deal('double');
            transientstreamlocs = table('Size',[length(allmaxlocs), length(transientstreamlocsvarnames)], 'VariableNames', transientstreamlocsvarnames, 'VariableTypes', transientstreamlocsvartypes);
            [transientstreamlocsvartypes{1:length(transientstreamlocsvarnames)}] = deal('double');
            transientstreamdata = zeros(length(allmaxlocs), (blstartsamples+posttransientsamples+1));
        end

        transientcount = 0;
        for eachmax = 1:length(allmaxlocs)
            % Initialize variables - set to empty
            currmaxloc = [];
            currmaxval = [];
            currpreminstartloc = [];
            currpreminendloc = [];   
            currminval = [];
            curramp = [];
            
            % Determine peak variables for inclusion
            currmaxloc = allmaxlocs(eachmax); % Index of current max in data stream
            currmaxval = data(eachfile).(whichstream)(currmaxloc); % Value of current max from data stream

            currpreminstartloc = currmaxloc-blstartsamples; % Find the start of the pre-peak baseline window
            currpreminendloc = currmaxloc-blendsamples; % Find the end of the pre-peak baseline window

            if currpreminstartloc < 1 % If the max is within the baseline length of the start of the session exclude it and move to the next max
                continue
            end

            currminval = min(data(eachfile).(whichstream)(currpreminstartloc:currpreminendloc)); % Find the value of pre-transient baseline window minimum
            currminloc = find(currminval == data(eachfile).(whichstream)(currpreminstartloc:currpreminendloc),1,'last'); % Index of the pre-transient baseline window min point in the data stream
            curramp = currmaxval-currminval; % Find the amplitude of the current transient relative to baseline

            if curramp >= data(eachfile).(whichthreshold) % If the current transient amplitude is greater than the threshold, quantify and add it to the table 'transientquantification'
                transientcount = transientcount + 1; % Update the total count of transients

                % Initialize variables - set to empty
                pretransientdata = [];
                curriseval = [];
                currrisesamples = [];
                currriseloc = [];
                posttransientdata = [];
                currfallval = [];
                currfallsamples = [];
                currfallloc = [];
                currAUCdata = [];
                currAUC = [];

                % Quantify rise time
                pretransientdata = data(eachfile).(whichstream)((currmaxloc-blstartsamples):currmaxloc); % Find pre-transient data
                currriseval = currminval + curramp*quantificationheight; % Find the rise value at the input quantification height (default is half height - 0.5)
                currrisesamples = length(pretransientdata) - find(pretransientdata >= currriseval,1,'first'); % Find the number of samples from rise start to transient peak
                currriseloc = currmaxloc-currrisesamples; % Find the location index of the rise start

                if (currmaxloc + posttransientsamples) > length(data(eachfile).(whichstream)) % Find post-transient data; check if post-transient period is within the length of the session
                    posttransientdata = data(eachfile).(whichstream)(currmaxloc:(currmaxloc+posttransientsamples));
                else
                    posttransientdata = data(eachfile).(whichstream)(currmaxloc:end); % If the end point of the post transient period is after the session end, just take to the end of the session
                end
                
                currfallval = currmaxval - curramp*quantificationheight; % Find the fall value at the input quantification height (default is half height - 0.5)
                currfallsamples = find(posttransientdata <= currfallval,1,'first'); % Find the number of samples from transient peak to fall end
                currfallloc = currmaxloc+currfallsamples; % Find the location index of the fall end
             
                if isempty(currfallsamples) % Catch for if no post-transient fall location is found
                    disp('WARNING: NO FALL FOUND FOR PEAK')
                    disp(transientcount)
                    disp('    Peak index: ')
                    disp(currmaxloc)
                    currfallval = NaN; % NaN out fall variables
                    currfallsamples = NaN;
                    currfallloc = NaN;
                    currwidthsamples = NaN;
                    currAUC = NaN;
                else % Calculate width and AUC (tranpezoidal method)
                    currwidthsamples = currfallloc - currriseloc;
                    currAUCdata = data(eachfile).(whichstream)(currriseloc:currfallloc);
                    currAUCdata = currAUCdata - min(currAUCdata);
                    currAUC = round(trapz(currAUCdata));
                end
                    
                % Add variables to the table 'transientquantification'
                transientquantification.transientID(transientcount) = transientcount;
                transientquantification.maxloc(transientcount) = currmaxloc;
                transientquantification.maxval(transientcount) = currmaxval;
                transientquantification.preminloc(transientcount) = currminloc;
                transientquantification.preminval(transientcount) = currminval;
                transientquantification.amp(transientcount) = curramp;
                transientquantification.risestartloc(transientcount) = currriseloc;
                transientquantification.risestartval(transientcount) = currriseval;
                transientquantification.risesamples(transientcount) = currrisesamples;
                transientquantification.risems(transientcount) = (currrisesamples/fs)*1000;            
                transientquantification.fallendloc(transientcount) = currfallloc;
                transientquantification.fallendval(transientcount) = currfallval;
                transientquantification.fallsamples(transientcount) = currfallsamples;
                transientquantification.fallms(transientcount) = (currfallsamples/fs)*1000;
                transientquantification.widthsamples(transientcount) = currwidthsamples;    
                transientquantification.widthms(transientcount) = (currwidthsamples/fs)*1000;
                transientquantification.AUC(transientcount) = currAUC;    

                if outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, add transient data stream
                    % Add locs for the cut data streams
                    transientstreamlocs.transientID(transientcount) = transientcount; 
                    transientstreamlocs.maxloc(transientcount) = blstartsamples;
                    transientstreamlocs.preminloc(transientcount) = blstartsamples-(currmaxloc-currminloc);
                    transientstreamlocs.risestartloc(transientcount) = blstartsamples-currrisesamples;
                    transientstreamlocs.fallendloc(transientcount) = blstartsamples+currfallsamples;
                    % Cut data streams
                    if (currmaxloc+posttransientsamples) <= length(data(eachfile).(whichstream))
                        transientstreamdata(transientcount,:) = data(eachfile).(whichstream)(currpreminstartloc:(currpreminstartloc+blstartsamples+posttransientsamples));
                    else
                        transientstreamdata(transientcount,:) = [data(eachfile).(whichstream)(currpreminstartloc:end),NaN(1,(blstartsamples+posttransientsamples)-length(data(eachfile).(whichstream)(currpreminstartloc:end)))];
                    end
                end
            end               
        end
        % Add inputs and transient quantification to the data structure
        data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).inputs = inputs;
        data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).transientquantification = transientquantification(1:transientcount,:);
        
        % OPTIONAL: If outputtransientdata is set to 1, add cut transient data streams and stream locs to data structure
        if outputtransientdata == 1
            data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).transientstreamlocs = transientstreamlocs(1:transientcount,:);
            data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).transientstreamdata = transientstreamdata(1:transientcount,:);
        end
    end
end
