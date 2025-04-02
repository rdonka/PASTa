function [data] = findSessionTransients_blmin(data,whichstream,whichthreshold,whichfs,preminstartms,preminendms,posttransientms,compoundtransientwindowms,quantificationheight,outputtransientdata)
% FINDSESSIONTRANSIENTS_BLMIN   Detects and quantifies transients in a data stream
%                               for the entire session using the pre-transient 
%                               baseline window minimum.
%
%   FINDSESSIONTRANSIENTS_BLMIN(DATA, WHICHBLTYPE, WHICHSTREAM, WHICHTHRESHOLD, WHICHFS...)
%   analyzes the specified data stream to detect transient events based on
%   the pre-transient baseline window minimum and specified threshold. This
%   function is called by the wrapper function FINDSESSIONTRANSIENTS.
%
% REQUIRED INPUTS:
%       DATA            - Structure array; each element corresponds to a session
%                         and must contain the fields specified by WHICHSTREAM,
%                         WHICHTHRESHOLD, and WHICHFS.
%
%       WHICHSTREAM     - String; name of the field in DATA containing the data stream
%                         to be analyzed (e.g., 'sigfiltz_normsession').
%
%       WHICHTHRESHOLD  - String; name of the field in DATA containing the numeric
%                         threshold values for transient detection (e.g., 'threshold_3SD').
%                         Thresholds should be precomputed and typically set to 2-3 standard deviations.
%
%       WHICHFS         - String; name of the field in DATA containing the sampling rate (fs)
%                         of the data stream.
%
%       PREMINSTARTMS             - Numeric; start time (ms) of the pre-transient baseline window.
%                                   Default: 800 ms.
%
%       PREMINDENDMS              - Numeric; end time (ms) of the pre-transient baseline window.
%                                   Default: 100 ms.
%
%       POSTTRANSIENTMS           - Numeric; duration (ms) after the transient peak for analysis.
%                                   Default: 2000 ms.
%
%       COMPOUNDTRANSIENTWINDOWMS - Numeric; window (ms) to search before and after each event
%                                   for compound transients. Default: 2000 ms.
%
%       QUANTIFICATIONHEIGHT      - Numeric; height (as a fraction of peak amplitude) at which to
%                                   characterize rise time, fall time, peak width, and area under
%                                   the curve (AUC). Must be between 0 and 1. Default: 0.5.
%
%       OUTPUTTRANSIENTDATA       - Logical; if true (1), outputs cut data streams for each transient
%                                   event. If false (0), skips this output. Default: true (1).
%
%   OUTPUTS:
%       DATA            - Structure array; each element corresponds to a session and includes
%                         the following added fields:
%                           - sessiontransients_blmin_<THRESHOLDLABEL>: A structure containing:
%                               - params: Structure of input parameters used for transient detection.
%                               - transientquantification: Table of quantified variables for each transient,
%                                 including amplitude, rise time, fall time, width, and AUC.
%                               - transientstreamlocs: Table of pre-transient baseline, transient peak,
%                                 rise, and fall locations for each transient.
%                               - transientstreamdata: Table of cut data streams from baseline start to
%                                 the end of the post-transient period for each transient event.
%
%
% See also: findSessionTransients, findSessionTransients_blmean, findSessionTransients_localmin
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Prepare Settings
% Import required and optional params into a structure
    params = struct(...
        'whichstream',whichstream,...
        'whichthreshold',whichthreshold,...
        'whichfs',whichfs,...
        'preminstartms', preminstartms,...
        'preminendms',preminendms,...
        'posttransientms',posttransientms,...
        'compoundtransientwindowms',compoundtransientwindowms,...
        'quantificationheight',quantificationheight,...
        'outputtransientdata',outputtransientdata);

    %% Find transients
    for eachfile = 1:length(data)
        disp(['Finding Transients: File ',num2str(eachfile)]) % Display which file is being processed
        try
            % Prep variables
            fs = data(eachfile).(whichfs);
            blstartsamples = floor(fs*(preminstartms/1000));
            blendsamples = floor(fs*(preminendms/1000));
            posttransientsamples = floor(fs*(posttransientms/1000));
            compoundtransientwindowsamples = floor(fs*(compoundtransientwindowms/1000));
            
            % FOR OUTPUT TRANSIENT DATA:
            premaxstart = floor(2*fs);
            postmaxend = floor(3*fs);
    
            % Find all maxes and mins in stream
            allmaxlocs = find(islocalmax(data(eachfile).(whichstream))); % Find all maxes
    
            % Prepare tables - preallocate size
            allvarnames = {'transientID','maxloc', 'maxval', 'preminloc', 'preminval', 'amp','risestartloc','risestartval','risesamples', 'risems',...
                'fallendloc','fallendval','fallsamples','fallms','widthsamples','widthms','AUC','IEIsamples','IEIs','compoundeventnum'};
            [allvartypes{1:length(allvarnames)}] = deal('double');
            transientquantification = table('Size',[length(allmaxlocs), length(allvarnames)], 'VariableNames', allvarnames, 'VariableTypes', allvartypes);
    
            if outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, prep data table and structure for individual transient streams
                transientstreamlocsvarnames = {'transientID','maxloc','preminloc','risestartloc','fallendloc'};
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
    
                currminval = min(data(eachfile).(whichstream)(currpreminstartloc:currpreminendloc)); % Find the value of pre-transient baseline window minimum
                currminloc = currpreminstartloc+find(currminval == data(eachfile).(whichstream)(currpreminstartloc:currpreminendloc),1,'last'); % Index of the pre-transient baseline window min point in the data stream
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
                    currIEIsamples = [];
                    currIEIs = [];
    
                    % Quantify rise time
                    pretransientdata = data(eachfile).(whichstream)((currmaxloc-blstartsamples):currmaxloc); % Find pre-transient data
                    currriseval = currminval + curramp*quantificationheight; % Find the rise value at the input quantification height (default is half height - 0.5)
                    currrisesamples = length(pretransientdata) - find(pretransientdata >= currriseval,1,'first'); % Find the number of samples from rise start to transient peak
                    currriseloc = currmaxloc-currrisesamples; % Find the location index of the rise start
    
                    if (currmaxloc + posttransientsamples) < length(data(eachfile).(whichstream)) % Find post-transient data; check if post-transient period is within the length of the session
                        posttransientdata = data(eachfile).(whichstream)(currmaxloc:(currmaxloc+posttransientsamples));
                    else
                        posttransientdata = data(eachfile).(whichstream)(currmaxloc:end); % If the end point of the post transient period is after the session end, just take to the end of the session
                    end
                    
                    % Quantify fall time and AUC
                    currfallval = currmaxval - curramp*quantificationheight; % Find the fall value at the input quantification height (default is half height - 0.5)
                    currfallsamples = find(posttransientdata <= currfallval,1,'first'); % Find the number of samples from transient peak to fall end
                    currfallloc = currmaxloc+currfallsamples; % Find the location index of the fall end
                 
                    if isempty(currfallsamples) % Catch for if no post-transient fall location is found
                        %disp('WARNING: NO FALL FOUND FOR PEAK')
                        %disp(transientcount)
                        %disp('    Peak index: ')
                        %disp(currmaxloc)
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
                        
                    % Find interevent interval (IEI)
                    if transientcount > 1
                        prevmaxloc = transientquantification.maxloc(transientcount-1);
                        currIEIsamples = currmaxloc - prevmaxloc;
                        currIEIs = currIEIsamples/fs;
                    else
                        currIEIsamples = NaN;
                        currIEIs = NaN;
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
                    transientquantification.IEIsamples(transientcount) = currIEIsamples;
                    transientquantification.IEIs(transientcount) = currIEIs;
    
                    if outputtransientdata == 1 % OPTIONAL: If outputtransientdata is set to 1, add transient data stream
                        % Add locs for the cut data streams
                        transientstreamlocs.transientID(transientcount) = transientcount; 
                        transientstreamlocs.maxloc(transientcount) = premaxstart;
                        transientstreamlocs.blstart(transientcount) = blstartsamples;
                        transientstreamlocs.blend(transientcount) = currmaxloc-blendsamples;
                        transientstreamlocs.preminloc(transientcount) = premaxstart-(currmaxloc-currminloc);
                        transientstreamlocs.risestartloc(transientcount) = premaxstart-currrisesamples;
                        transientstreamlocs.fallendloc(transientcount) = premaxstart+currfallsamples;
                        
                        % Cut data streams
                        if (currmaxloc-premaxstart) >= 1 && (currmaxloc+postmaxend) <= length(data(eachfile).(whichstream))
                            transientstreamdata(transientcount,:) = data(eachfile).(whichstream)(currmaxloc-premaxstart:currmaxloc+postmaxend);
                        else
                            if (currmaxloc-premaxstart) < 1
                                premissingsamples = abs(currmaxloc-premaxstart);
                                prenans = NaN(1,premissingsamples);
                            else
                                premissingsamples = 0;
                                prenans = [];
                            end
                            if (currmaxloc+postmaxend) >= length(data(eachfile).(whichstream))
                                postmissingsamples = abs((postmaxend)-length(data(eachfile).(whichstream)(currmaxloc:end)));
                                postnans = [NaN(1,postmissingsamples)];
                            else
                                postmissingsamples = 0;
                                postnans = [];
                            end
                            streamwithnans = [prenans,data(eachfile).(whichstream)(currmaxloc-(premaxstart-premissingsamples-1):(postmaxend-postmissingsamples)),postnans];
                            streamwithnans = [streamwithnans, NaN(1,(premaxstart+postmaxend+1)-length(streamwithnans))];
                            transientstreamdata(transientcount,:) = streamwithnans;
                        end
                    end
                end               
            end

            % Check for compound peak. 0 = no compound peaks in window; Otherwise, value reflects the event number of the transient relative to the others in the window.
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
            data(eachfile).(append('sessiontransients_blmin_',whichthreshold)).params = params;
            data(eachfile).(append('sessiontransients_blmin_',whichthreshold)).transientquantification = transientquantification(1:transientcount,:);

            % Display how many transients were found
            disp(append('   Total Transients: ',num2str(transientcount)))
            
            % OPTIONAL: If outputtransientdata is set to 1, add cut transient data streams and stream locs to data structure
            if outputtransientdata == 1
                data(eachfile).(append('sessiontransients_blmin_',whichthreshold)).transientstreamlocs = transientstreamlocs(1:transientcount,:);
                data(eachfile).(append('sessiontransients_blmin_',whichthreshold)).transientstreamdata = transientstreamdata(1:transientcount,:);
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