function [data] = findSessionTransients_blmean(data,whichstream,whichthreshold,whichfs,thresholdlabel,varargin)
% FINDSESSIONTRANSIENTS   Finds transients for the whole session of data.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the data
%                       stream you want to analyze.
%
%       WHICHSTREAM:    A variable containing a string with the name of the 
%                       field containing the stream to be analyzed for 
%                       transients.
%
%       WHICHTHRESHOLD: A variable containing a string with the name of the 
%                       field containing the threshold criterion for the 
%                       stream to be analyzed for transients.
%
%       WHICHBLSTART:   The name of the field containing the number of
%                       samples before the pk to include in the baseline 
%                       window.
%
%       WHICHBLEND:     The name of the field containing the number of
%                       samples before the pk to exclude from the baseline 
%                       window.
%
%       WHICHFS:        The name of the field containing the sampling rate
%                       of the streams (fs).
%
% OPTIONAL INPUTS:
%       QUANTIFICATIONHEIGHT:         The height at which to characterize rise time, fall
%                       time, and AUC. Must be a number between 0 and 1.
%                       Default: 0.5
% OUTPUTS:
%       DATA:           The original data structure with
%                       sessiontransients_blmean_THRESHOLDLABEL added in.
%                       Output includes the max peak location and value, 
%                       baseline min start and end location, baseline 
%                       window mean, and peak amplitude. Peak rise and fall 
%                       time are determined at the specified height 
%                       (default is half height) and output in samples. 
%                       AUC is determined at the specified height.
%
% Written by R M Donka, October 2024.
% Stored in PASTa GitHub repository, see Wiki for additional notes.

disp('SESSION TRANSIENTS: Peak baseline determined by mean value of the specified baseline window.')

%% Prepare Settings
% Import optional inputs into a structure
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

    if isempty(inputs.preminstartms)
        preminstartms = 1000;
        inputs.preminstartms = preminstartms;
    else
        preminstartms = inputs.preminstartms;
    end

    if isempty(inputs.preminendms)
        preminendms = 100;
        inputs.preminendms = preminendms;
    else
        preminendms = inputs.preminendms;
    end
    if isempty(inputs.quantificationheight)
        quantificationheight = 0.5; % Height for rise/fall/AUC defaults to 0.5 (half height)
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
    if isempty(inputs.posttransientms)
        posttransientms = 2000; % ms post transient max for fall/AUC and output data
        inputs.posttransientms = posttransientms;
    else
        posttransientms = inputs.posttransientms;
    end

    disp('INPUTS:') % Display all input values
    disp(inputs)
        
    for eachfile = 1:length(data)
        disp(['Finding Transients: File ',num2str(eachfile)]) % Display which file is being processed
        
        % Prep variables
        fs = data(eachfile).(whichfs);
        blstartsamples = floor(fs*(preminstartms/1000));
        blendsamples = floor(fs*(preminendms/1000));
        posttransientsamples = floor(fs*(posttransientms/1000));

        % Find maxes in stream
        allmaxlocs = find(islocalmax(data(eachfile).(whichstream))); % Find all maxes
        
        % Prepare tables - preallocate size
        allvarnames = {'transientID','maxloc', 'maxval', 'preminstartloc', 'preminendloc', 'preminval', 'amp', 'risestartloc','risesamples', 'risems','fallendloc', 'fallsamples', 'fallms', 'AUC'};
        [allvartypes{1:length(allvarnames)}] = deal('double');
        transientquantification = table('Size',[length(allmaxlocs), length(allvarnames)], 'VariableNames', allvarnames, 'VariableTypes', allvartypes);

        if outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, prep data table and structure for individual transient streams
            transientstreamlocsvarnames = {'transientID','maxloc','currpreminstartloc', 'currpreminendloc','posttransientloc','risestartloc','fallendloc'};
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
            
            currmaxloc = allmaxlocs(eachmax); % Index of current max in data stream
            currmaxval = data(eachfile).(whichstream)(currmaxloc); % Value of current max from data stream

            currpreminstartloc = currmaxloc-blstartsamples; % Find the start of the pre-peak baseline window
            currpreminendloc = currmaxloc-blendsamples; % Find the end of the pre-peak baseline window

            if currpreminstartloc < 1
                continue
            end

            currpreminstartloc = currmaxloc-blstartsamples; % Find the start of the pre-peak baseline window
            currpreminendloc = currmaxloc-blendsamples; % Find the end of the pre-peak baseline window

            currminval = mean(data(eachfile).(whichstream)(currpreminstartloc:currpreminendloc)); % Find mean of baseline window
            curramp = currmaxval-currminval; % Current peak amplitude

            if curramp >= data(eachfile).(whichthreshold) % If the current peak is greater than the threshold, add to the table transientquantification
                transientcount = transientcount + 1;

                % Initialize variables - set to empty
                prepkdata = [];
                curriseval = [];
                currrisesamples = [];
                currriseloc = [];
                postpkdata = [];
                currfallval = [];
                currfallsamples = [];
                currfallloc = [];
                currpkAUCdata = [];
                currAUC = [];

                prepkdata = data(eachfile).(whichstream)((currmaxloc-blstartsamples):currmaxloc);
                currriseval = currminval + curramp*quantificationheight;
                currrisesamples = length(prepkdata) - find(prepkdata >= currriseval,1,'first');
                currriseloc = currmaxloc-currrisesamples;

                if (currmaxloc + posttransientsamples) > length(data(eachfile).(whichstream))
                    postpkdata = data(eachfile).(whichstream)(currmaxloc:(currmaxloc+posttransientsamples));
                else
                    postpkdata = data(eachfile).(whichstream)(currmaxloc:end);
                end
                
                currfallval = currmaxval - curramp*quantificationheight;
                currfallsamples = find(postpkdata <= currfallval,1,'first');
                currfallloc = currmaxloc+currfallsamples;
                
                if isempty(currfallsamples)
                    disp('WARNING: NO FALL FOUND FOR PEAK')
                    disp(transientcount)
                    disp('    Peak index: ')
                    disp(currmaxloc)
                    currfallval = NaN;
                    currfallsamples = NaN;
                    currfallloc = NaN;
                    currAUC = NaN;
                else
                    currpkAUCdata = data(eachfile).(whichstream)(currriseloc:currfallloc);
                    currpkAUCdata = currpkAUCdata - min(currpkAUCdata);
                    currAUC = round(trapz(currpkAUCdata));
                end
                    
                transientquantification.transientID(transientcount) = transientcount;
                transientquantification.maxloc(transientcount) = currmaxloc;
                transientquantification.maxval(transientcount) = currmaxval;
                transientquantification.preminstartloc(transientcount) = currpreminstartloc;
                transientquantification.preminendloc(transientcount) = currpreminendloc;
                transientquantification.preminval(transientcount) = currminval;
                transientquantification.amp(transientcount) = curramp;
                transientquantification.risestartloc(transientcount) = currriseloc;
                transientquantification.risesamples(transientcount) = currrisesamples;
                transientquantification.risems(transientcount) = (currrisesamples/fs)*1000;            
                transientquantification.fallendloc(transientcount) = currfallloc;
                transientquantification.fallsamples(transientcount) = currfallsamples;
                transientquantification.fallms(transientcount) = (currfallsamples/fs)*1000;                
                transientquantification.AUC(transientcount) = currAUC;    


                if outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, add transient data stream
                    transientstreamlocs.transientID(transientcount) = transientcount;
                    transientstreamlocs.maxloc(transientcount) = blstartsamples;
                    transientstreamlocs.preminstartloc(transientcount) = 1;
                    transientstreamlocs.preminendloc(transientcount) = blstartsamples-blendsamples;
                    transientstreamlocs.posttransientloc(transientcount) = blstartsamples+posttransientsamples;
                    transientstreamlocs.risestartloc(transientcount) = blstartsamples-currrisesamples;
                    transientstreamlocs.fallendloc(transientcount) = blstartsamples+currfallsamples;

                    if (currmaxloc+posttransientsamples) <= length(data(eachfile).(whichstream))
                        transientstreamdata(transientcount,:) = data(eachfile).(whichstream)(currpreminstartloc:(currpreminstartloc+blstartsamples+posttransientsamples));
                    else
                        transientstreamdata(transientcount,:) = [data(eachfile).(whichstream)(currpreminstartloc:end),NaN(1,(blstartsamples+posttransientsamples)-length(data(eachfile).(whichstream)(currpreminstartloc:end)))];
                    end
                end
            end               
        end

        data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).inputs = inputs;
        data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).transientquantification = transientquantification(1:transientcount,:);
        
        if outputtransientdata == 1
            data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).transientstreamlocs = transientstreamlocs(1:transientcount,:);
            data(eachfile).(append('sessiontransients_blmin_',thresholdlabel)).transientstreamdata = transientstreamdata(1:transientcount,:);
        end
    end
end

