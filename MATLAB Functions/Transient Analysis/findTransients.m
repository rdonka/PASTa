function [transientdata] = findTransients(data,addvariablesfieldnames,streamfieldname,thresholdfieldname,fsfieldname,varargin)
% FINDTRANSIENTS     Detects and quantifies transients in a data stream
%                           for the entire session using the pre-transient 
%                           baseline window mean.
%
%   FINDTRANSIENTS(DATA,STREAMFIELDNAME, THRESHOLDFIELDNAME, FSFIELDNAME, 'PARAM1', VAL1, ...)
%   analyzes the specified data stream to detect transient events based on
%   the specified threshold.
%
% REQUIRED INPUTS:
%       DATA                - Structure array; each element corresponds to a session
%                             and must contain the fields specified by STREAMFIELDNAME,
%                             THRESHOLDFIELDNAME, and FSFIELDNAME.
%
%       ADDVARIABLESFIELDNAMES - Cell array; names of the fields in DATA to
%                             add to the new data structure. This should include SubjectID 
%                             and experimentally relevant metadata.
%                             (e.g., {'SubjectID', 'BlockFolder', 'Dose'})
%
%       STREAMFIELDNAME     - String; name of the field in DATA containing the data stream
%                             to be analyzed (e.g., 'sigfiltz_normsession').
%
%       THRESHOLDFIELDNAME  - String; name of the field in DATA containing the numeric
%                             threshold values for transient detection (e.g., 'SDthreshold').
%                             Thresholds should be precomputed and typically 
%                             set to 2.6 standard deviations.
%
%       FSFIELDNAME         - String; name of the field in DATA containing the sampling rate (fs)
%                             of the data stream.
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%       'bltype'                  - String; Method for pre-transient peak baseline detection. 
%                                   Options are 'blmean', 'blmin', and 'localmin'.
%                                   Default: 'blmean'.
%
%       'blstartms'           - Numeric; start time (ms) of the pre-transient baseline window.
%                                   Default: 1000 ms.
%
%       'blendms'             - Numeric; end time (ms) of the pre-transient baseline window.
%                                   Default: 200 ms.
%
%       'posttransientms'         - Numeric; duration (ms) after the transient peak for analysis.
%                                   Default: 2000 ms.
%
%
%       'quantificationheight'    - Numeric; height (as a fraction of peak amplitude) at which to
%                                   characterize rise time, fall time, peak width, and area under
%                                   the curve (AUC). Must be between 0 and 1. Default: 0.5.
%
%       'compoundtransientwindowms' - Numeric; window (ms) to search before and after each event
%                                     for compound transients. Default: 2000 ms.
%
%       'AUCwindowms'         - Numeric; window (ms) for transient AUC
%                               window.
%
%       'outputtransientdata'     - Logical; if true (1), outputs cut data streams for each transient
%                                   event. If false (0), skips this output.
%                                   Default: true (1).
%
%       'outputpremaxS'           - Numeric; Number of seconds pre transient peak to include in 
%                                   transient data output streams. Default: 3.
%
%       'outputpostmaxS'           - Numeric; Number of seconds post transient peak to include in 
%                                   transient data output streams. Default: 5.
%
%   OUTPUTS:
%       TRANSIENTDATA    - Structure array; each element corresponds to a session and includes 
%                          all fields specified in ADDVARIABLESFIELDNAMES as well as the following:
%                           - params.findTransients: Structure of input parameters used for transient detection.
%                           - transientquantification: Table of quantified variables for each transient,
%                             including amplitude, rise time, fall time, width, and AUC.
%                           - If OUTPUTTRANSIENTDATA is set to 1:
%                               - transientstreamlocs: Table of pre-transient baseline, transient peak,
%                                 rise, and fall locations for each transient.
%                               - transientstreamdata: Table of cut data streams from baseline start to
%                                 the end of the post-transient period for each transient event.
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Prepare Settings
    % Prepare default values
    defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.

    % Import required and optional params into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details

    % Add optional name-value pair arguments with validation
    addParameter(p, 'bltype', defaultparameters.bltype, @(x) ischar(x) && ismember(x, {'blmean', 'blmin', 'localmin'})); % bltype: Pre-peak baseline method
    addParameter(p, 'blstartms', defaultparameters.blstartms, @isnumeric); % blstartms: Numeric (ms); pre-transient to start the baseline period
    addParameter(p, 'blendms', defaultparameters.blendms, @isnumeric); % blendms: Numeric (ms); ms pre-transient to start the baseline period
    addParameter(p, 'posttransientms', defaultparameters.posttransientms, @isnumeric); % posttransientms: Numeric (ms); ms post-transient to use as the oost-transient fall window
    addParameter(p, 'AUCwindowms',  defaultparameters.compoundtransientwindowms, @isnumeric);  % Numeric (ms); Window size to search before and after each event for compound transients
    addParameter(p, 'compoundtransientwindowms',  defaultparameters.compoundtransientwindowms, @isnumeric);  % Numeric (ms); Window size to search before and after each event for compound transients
    addParameter(p, 'quantificationheight', defaultparameters.quantificationheight, @(x) isnumeric(x) && x >= 0 && x <= 1); % quantificationheight: Numeric between 0 and 1; Height for quantification of transients (rise/fall/AUC)
    addParameter(p, 'outputtransientdata', defaultparameters.outputtransientdata, @islogical); % outputtransientdata: Logical; If set to 1 (true), data streams for individual transients will be added to the data structure
    addParameter(p, 'outputpremaxS', defaultparameters.outputpremaxS, @isnumeric); % outputpremaxS: Numeric; Number of seconds pre transient peak to include in transient data output
    addParameter(p, 'outputpostmaxS', defaultparameters.outputpostmaxS, @isnumeric); % outputpostmaxS: Logical; Numeric; Number of seconds post transient peak to include in transient data output

    parse(p, varargin{:});

    % Retrieve parsed params into params structure
    params = p.Results;

    % Add function call variables
    params.streamfield = streamfieldname;
    params.thresholdfield = thresholdfieldname;
    params.fsfield = fsfieldname;

    %% Display settings
    switch params.bltype
        case 'blmean'
            disp(['FINDTRANSIENTS: Identifying and quantifying transients in stream ',streamfieldname,' with thresholds from ',thresholdfieldname,'. Pre-peak baselines will be determined by the baseline window mean.'])
        case 'blmin'
            disp(['FINDTRANSIENTS: Identifying and quantifying transients in stream ',streamfieldname,' with thresholds from ',thresholdfieldname,'. Pre-peak baselines will be determined by the baseline window minimum.'])
        case 'localmin'
            disp(['FINDTRANSIENTS: Identifying and quantifying transients in stream ',streamfieldname,' with thresholds from ',thresholdfieldname,'. Pre-peak baselines will be determined by the local minimum.'])
        otherwise
            error('No viable transient baseline type specific. BLTYPE must be set to blmin, blmean, or localmin.')
    end
    disp('   Data will be output to a new data structure with the following variables from addvariablesfieldnames:')
    disp(addvariablesfieldnames)
    disp('   PARAMETERS:') % Display all input values
    disp(params)

    % Prepare transient output structure
    transientdata = struct();
    for eachvariable = 1:length(addvariablesfieldnames)
        currvariablefieldname = char(addvariablesfieldnames(eachvariable));
        for eachfile = 1:length(data)
            transientdata(eachfile).(currvariablefieldname) = [data(eachfile).(currvariablefieldname)];
        end
    end


    %% Find transients
    for eachfile = 1:length(data)
        disp(['Finding Transients: File ',num2str(eachfile)]) % Display which file is being processed
        try
            % Prep variables
            datastream = data(eachfile).(streamfieldname);
            fs = data(eachfile).(fsfieldname);
            blstartsamples = floor(fs*(params.blstartms/1000));
            blendsamples = floor(fs*(params.blendms/1000));
            posttransientsamples = ceil(fs*(params.posttransientms/1000));
            compoundtransientwindowsamples = floor(fs*(params.compoundtransientwindowms/1000));
            
            % Find all maxes in stream
            allmaxlocs = find(islocalmax(datastream)); % Find all maxes
            
            % Prepare tables - preallocate size
            allvarnames = {'transientID','maxloc', 'maxval', 'blstartloc', 'blendloc', 'blloc', 'blval', 'amp','quantheightval','risestartloc','risesamples', 'risems',...
                'fallendloc','fallsamples','fallms','widthsamples','widthms','AUC','IEIsamples','IEIms','IEIs','compoundeventnum'};
            [allvartypes{1:length(allvarnames)}] = deal('double');
            transientquantification = table('Size',[length(allmaxlocs), length(allvarnames)], 'VariableNames', allvarnames, 'VariableTypes', allvartypes);
    
            if params.outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, prep data table and structure for individual transient streams
                % Prepare indexes to cut streams
                premaxstart = floor(params.outputpremaxS*fs);
                postmaxend = ceil(params.outputpostmaxS*fs);

                transientstreamlocsvarnames = {'transientID','maxloc','blstartloc', 'blendloc','risestartloc','fallendloc'};
                [transientstreamlocsvartypes{1:length(transientstreamlocsvarnames)}] = deal('double');
                transientstreamlocs = table('Size',[length(allmaxlocs), length(transientstreamlocsvarnames)], 'VariableNames', transientstreamlocsvarnames, 'VariableTypes', transientstreamlocsvartypes);
                [transientstreamlocsvartypes{1:length(transientstreamlocsvarnames)}] = deal('double');
                transientstreamdata = zeros(length(allmaxlocs), (premaxstart+postmaxend+1));
            end
    
            transientcount = 0;

            for eachmax = 1:length(allmaxlocs)
                % Initialize variables - set to empty
                currmaxloc = [];
                currmaxval = [];
                currblstartloc = [];
                currblendloc = [];   
                currminval = [];
                currminloc = [];
                curramp = [];
                
                % Determine peak variables for inclusion
                currmaxloc = allmaxlocs(eachmax); % Index of current max in data stream
                currmaxval = datastream(currmaxloc); % Value of current max from data stream

                % Find pre-peak baseline window indexes
                currblstartloc = currmaxloc-blstartsamples; % Find the start of the pre-peak baseline window
                currblendloc = currmaxloc-blendsamples; % Find the end of the pre-peak baseline window
    
                if currblstartloc < 1 % If the max is within the baseline length of the start of the session exclude it and move to the next max
                    continue
                end

                % Find pre-peak baseline
                switch params.bltype
                    case 'blmean'
                        currminloc = currblendloc - ((blstartsamples-blendsamples)/2);
                        currminval = mean(datastream(currblstartloc:currblendloc)); % Find mean of pre-transient baseline window
                    case 'blmin'
                        currminval = min(datastream(currblstartloc:currblendloc)); % Find the value of pre-transient baseline window minimum
                        currminloc = currblstartloc+find(currminval == datastream(currblstartloc:currblendloc),1,'last'); % Index of the pre-transient baseline window min point in the data stream
                    case 'localmin'
                        currminloc = find(islocalmin(datastream(currblstartloc:currblendloc)),1,'last') + currblstartloc; % Find the last local minimum in the baseline window
                        currminval = datastream(currminloc); % Find the value of the pre-transient local minimum
                end

                curramp = currmaxval-currminval; % Find the amplitude of the current transient relative to baseline
    
                if curramp >= data(eachfile).(thresholdfieldname) % If the current transient amplitude is greater than the threshold, quantify and add it to the table 'transientquantification'
                    transientcount = transientcount + 1; % Update the total count of transients
   
                    % Initialize variables - set to empty
                    pretransientdata = [];
                    quantheightval = [];
                    currrisesamples = [];
                    currriseloc = [];
                    posttransientdata = [];
                    currfallsamples = [];
                    currfallloc = [];
                    currpkAUCwindowdata = [];
                    currAUCwindow = [];
                    currpkAUCdata = [];
                    currAUC = [];
                    currIEIsamples = [];
                    currIEIms = [];
                    currIEIs = [];
    

                    % Quantify rise time
                    pretransientdata = datastream((currmaxloc-blstartsamples):currmaxloc); % Extract pre-transient data from baseline start to peak max
                    quantheightval = currmaxval - curramp*params.quantificationheight; % Find the rise value at the input quantification height
                    currrisesamples = length(pretransientdata) - find(pretransientdata >= quantheightval,1,'first'); % Find the number of samples from rise start to transient peak
                    currrisems = (currrisesamples/fs)*1000;
                    currriseloc = currmaxloc-currrisesamples; % Find the location index of the rise start
    
                    % Quantify fall time
                    if (currmaxloc + posttransientsamples) < length(datastream) % Find post-transient data; check if post-transient period is within the length of the session
                        posttransientdata = datastream(currmaxloc:(currmaxloc+posttransientsamples));
                    else
                        posttransientdata = datastream(currmaxloc:end); % If the end point of the post transient period is after the session end, just take to the end of the session
                    end
                    
                    currfallsamples = find(posttransientdata <= quantheightval,1,'first'); % Find the number of samples from transient peak to fall end
                    currfallloc = currmaxloc+currfallsamples; % Find the location index of the fall end
                    
                    %  Quantify transient AUC - window
                    AUCwindowsamples = floor((params.AUCwindowms/1000)*fs);
                    currpkAUCwindowstart = currriseloc - (currriseloc - currblendloc);
                    currpkAUCwindowend = currriseloc+AUCwindowsamples;
                    if currpkAUCwindowend > length(datastream)
                        currAUCwindow = NaN;
                    else
                        currpkAUCwindowdata = [datastream(currpkAUCwindowstart:currpkAUCwindowend)] - min(datastream(currpkAUCwindowstart:currpkAUCwindowend));
                        currAUCwindow = round(trapz(currpkAUCwindowdata));
                    end

                    if isempty(currfallsamples) % Catch for if no post-transient fall location is found
                        currfallsamples = NaN;
                        currfallloc = NaN;
                        currfallms = NaN;
                        currwidthsamples = NaN;
                        currwidthms = NaN;
                        currAUC = NaN;
                    else % Calculate width and AUC (tranpezoidal method)
                        currfallms = (currfallsamples/fs)*1000;
                        currwidthsamples = currfallloc - currriseloc;
                        currwidthms =(currwidthsamples/fs)*1000;
                        currpkAUCdata = [datastream(currriseloc:currfallloc)] - min(datastream(currriseloc:currfallloc));
                        currAUC = round(trapz(currpkAUCdata));
                    end
                        
                    % Find interevent interval (IEI)
                    if transientcount > 1
                        prevmaxloc = transientquantification.maxloc(transientcount-1);
                        currIEIsamples = currmaxloc - prevmaxloc;
                        currIEIms = (currIEIsamples/fs)*1000;
                        currIEIs = currIEIsamples/fs;
                    else
                        currIEIsamples = NaN;
                        currIEIms = NaN;
                        currIEIs = NaN;
                    end

                    % Add variables to the table 'transientquantification'
                    transientquantification.transientID(transientcount) = transientcount;
                    transientquantification.maxloc(transientcount) = currmaxloc;
                    transientquantification.maxval(transientcount) = currmaxval;
                    transientquantification.blstartloc(transientcount) = currblstartloc;
                    transientquantification.blendloc(transientcount) = currblendloc;
                    transientquantification.blloc(transientcount) = currminloc;
                    transientquantification.blval(transientcount) = currminval;
                    transientquantification.amp(transientcount) = curramp;
                    transientquantification.quantheightval(transientcount) = quantheightval;
                    transientquantification.risestartloc(transientcount) = currriseloc;
                    transientquantification.risesamples(transientcount) = currrisesamples;
                    transientquantification.risems(transientcount) = currrisems;            
                    transientquantification.fallendloc(transientcount) = currfallloc;
                    transientquantification.fallsamples(transientcount) = currfallsamples;
                    transientquantification.fallms(transientcount) = currfallms;
                    transientquantification.widthsamples(transientcount) = currwidthsamples;    
                    transientquantification.widthms(transientcount) = currwidthms;
                    transientquantification.AUC(transientcount) = currAUC;
                    transientquantification.AUCwindow(transientcount) = currAUCwindow;
                    transientquantification.IEIsamples(transientcount) = currIEIsamples;
                    transientquantification.IEIms(transientcount) = currIEIms;
                    transientquantification.IEIs(transientcount) = currIEIs;
    
                    if params.outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, add transient data stream
                        % Add locs for the cut data streams
                        transientstreamlocs.transientID(transientcount) = transientcount; 
                        transientstreamlocs.maxloc(transientcount) = premaxstart;
                        transientstreamlocs.blstart(transientcount) = blstartsamples;
                        transientstreamlocs.blend(transientcount) = currmaxloc-blendsamples;
                        transientstreamlocs.risestartloc(transientcount) = premaxstart-currrisesamples;
                        transientstreamlocs.fallendloc(transientcount) = premaxstart+currfallsamples;
                        
                        % Cut data streams
                        if (currmaxloc-premaxstart) >= 1 && (currmaxloc+postmaxend) <= length(datastream)
                            transientstreamdata(transientcount,:) = datastream(currmaxloc-premaxstart:currmaxloc+postmaxend);
                        else
                            if (currmaxloc-premaxstart) < 1
                                premissingsamples = abs(currmaxloc-premaxstart);
                                prenans = NaN(1,premissingsamples);
                            else
                                premissingsamples = 0;
                                prenans = [];
                            end
                            if (currmaxloc+postmaxend) >= length(datastream)
                                postmissingsamples = abs((postmaxend)-length(datastream(currmaxloc:end)));
                                postnans = [NaN(1,postmissingsamples)];
                            else
                                postmissingsamples = 0;
                                postnans = [];
                            end
                            streamwithnans = [prenans,datastream(currmaxloc-(premaxstart-premissingsamples-1):(postmaxend-postmissingsamples-1)),postnans];
                            streamwithnans = [streamwithnans, NaN(1,(premaxstart+postmaxend+1)-length(streamwithnans))];
                            transientstreamdata(transientcount,:) = streamwithnans;
                        end
                    end
                end               
            end
            
            % Identify compound transients. 0 = no compound peaks in window; Otherwise, value reflects the event number of the transient relative to the others in the window.
            for eachtransient = 1:transientcount
                compoundstart = transientquantification{eachtransient,"maxloc"}-compoundtransientwindowsamples;
                compoundend = transientquantification{eachtransient,"maxloc"}+compoundtransientwindowsamples;

                transientsbefore = sum(transientquantification{1:transientcount,"maxloc"}>compoundstart & transientquantification{1:transientcount,"maxloc"}<transientquantification{eachtransient,"maxloc"});
                transientsafter = sum(transientquantification{1:transientcount,"maxloc"}>transientquantification{eachtransient,"maxloc"} & transientquantification{1:transientcount,"maxloc"}<compoundend);

                if (transientsbefore == 0) && (transientsafter == 0)
                    currcompoundeventnum = 0;
                else 
                    currcompoundeventnum = transientsbefore+1;
                end
                transientquantification.compoundeventnum(eachtransient) = currcompoundeventnum;
            end

            % Add params and transient quantification to the data structure
            transientdata(eachfile).params.(mfilename) = params; % Add overall parameters
            transientdata(eachfile).params.(mfilename).fs = fs;
            transientdata(eachfile).params.(mfilename).streamtotalsamples = length(datastream);
            transientdata(eachfile).transientquantification = transientquantification(1:transientcount,:);
            
            % Display how many transients were found
            disp(append('   Total Transients: ',num2str(transientcount)))
            
            % OPTIONAL: If outputtransientdata is set to 1, add cut transient data streams and stream locs to data structure
            if params.outputtransientdata == 1
                transientdata(eachfile).transientstreamlocs = transientstreamlocs(1:transientcount,:);
                transientdata(eachfile).transientstreamdata = transientstreamdata(1:transientcount,:);
            end
        catch ME
            fprintf('   ERROR: %s\n', ME.message);
        end
    end
end

% Copyright (C) 2025 Rachel Donka
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.