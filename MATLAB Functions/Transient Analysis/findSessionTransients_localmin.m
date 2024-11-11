function [data] = findSessionTransients_localmin(data,whichstream,whichthreshold,whichfs,preminstartms,preminendms,posttransientms,quantificationheight,outputtransientdata)
% FINDSESSIONTRANSIENTS_LOCALMIN  Finds transients for the whole session. 
%                               Pre-transient baselines are set to the
%                               local minimum directly preceding the
%                               transient within the baseline window. 
%
%                               NOTE: This sub-function is called by the
%                               main function FINDSESSIONTRANSIENTS. If
%                               called outside the main function, users
%                               must specify all input values manually.
%
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the data
%                       stream you want to analyze.
%
%       WHICHSTREAM:    A variable containing a string with the name of the 
%                       field containing the stream to be analyzed for 
%                       transients. For example, 'sigfiltz_normsession'.
%
%       WHICHTHRESHOLD: A variable containing a string with the name of the
%                       field containing the prepared numeric threshold
%                       values for each stream. For example, 'threshold_3SD'.
%                       NOTE: Threshold values should be calculated
%                       before using the findSessionTransients functions.
%                       Typically thresholds are set to 2-3 SDs. If the
%                       input data stream is Z scored, this can be the
%                       actual SD threshold number. If the input data
%                       stream is not Z scored, find the corresponding
%                       value to 2-3 SDs for each subject.
%                     
%       WHICHFS:        The name of the field containing the sampling rate
%                       of the streams (fs).
%
%       PREMINSTARTMS:  Number of millseconds pre-transient to use as the
%                       start of the baseline window. 
%                       NOTE: Default is set to 800 in the main
%                       'findSessionTransients' function.
%
%       PREMINENDMS:    Number of millseconds pre-transient to use as the
%                       end of the baseline window. If you want to
%                       determine the absolute local minimum before the
%                       transient, set this to 0.
%                       NOTE: Default is set to 100 in the main
%                       'findSessionTransients' function.
%
%       POSTTRANSIENTMS: Number of millseconds post-transient to use for
%                       the post peak baseline and trimmed data output.
%                       NOTE: Default is set to 2000 in the main
%                       'findSessionTransients' function.
%
%       QUANTIFICATIONHEIGHT: The height at which to characterize rise time,
%                       fall time, and AUC. Must be a number between 0 and 1.
%                       NOTE: Default is set to 0.5 in the main
%                       'findSessionTransients' function.
%  
%       OUTPUTTRANSIENTDATA: Set to 1 to output cut data streams for each
%                       transient event. Set to 0 to skip.
%                       NOTE: Default is set to 1 in the main
%                       'findSessionTransients' function.
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


%% Prepare Settings
% Import required and optional inputs into a structure
    inputs = struct(...
        'whichstream',whichstream,...
        'whichthreshold',whichthreshold,...
        'whichfs',whichfs,...
        'preminstartms', preminstartms,...
        'preminendms',preminendms,...
        'posttransientms',posttransientms,...
        'quantificationheight',quantificationheight,...
        'outputtransientdata',outputtransientdata);

    % Display settings
    disp("FIND SESSION TRANSIENTS: Peak baseline determined by last local minimum in the specified baseline window. WHICHBLTYPE set to 'localmin'")
    disp("     SUBFUNCTION: findSessionTransients_localmin")
    disp(append("     Transient data will be added to data structure as 'sessiontransients_localmin_",whichthreshold,"'." ))

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
        allminlocs = find(islocalmin(data(eachfile).(whichstream))); % Find all mins

        if allmaxlocs(1) < allminlocs(1) % If the first max is without a min, remove it
            allmaxlocs = allmaxlocs(2:end);
        end

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
            currminloc = [];
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

            allpreminlocs = allminlocs(allminlocs > currpreminstartloc & allminlocs < currpreminendloc); % Indexes of all mins that are between the baseline start and end point relative to the current max

            currminloc = allpreminlocs(find(allpreminlocs < allmaxlocs(eachmax),1,'last')); % Index of the pre-transient local minimum (last min point in the baseline window) in the data stream
            currminval = data(eachfile).(whichstream)(currminloc); % Find the value of pre-transient local minimum
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
                currpkAUCdata = [];
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
                    currpkAUCdata = data(eachfile).(whichstream)(currriseloc:currfallloc);
                    currpkAUCdata = currpkAUCdata - min(currpkAUCdata);
                    currAUC = round(trapz(currpkAUCdata));
                end
                    
                % Add variables to the table 'transientquantification'
                transientquantification.transientID(transientcount) = transientcount;
                transientquantification.maxloc(transientcount) = currmaxloc;
                transientquantification.maxval(transientcount) = currmaxval;
                transientquantification.preminstartloc(transientcount) = currpreminstartloc;
                transientquantification.preminendloc(transientcount) = currpreminendloc;
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
        data(eachfile).(append('sessiontransients_localmin_',whichthreshold)).inputs = inputs;
        data(eachfile).(append('sessiontransients_localmin_',whichthreshold)).transientquantification = transientquantification(1:transientcount,:);
        
        % OPTIONAL: If outputtransientdata is set to 1, add cut transient data streams and stream locs to data structure
        if outputtransientdata == 1
            data(eachfile).(append('sessiontransients_localmin_',whichthreshold)).transientstreamlocs = transientstreamlocs(1:transientcount,:);
            data(eachfile).(append('sessiontransients_localmin_',whichthreshold)).transientstreamdata = transientstreamdata(1:transientcount,:);
        end
    end
end

% Copyright (C) 2024 Rachel Donka
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