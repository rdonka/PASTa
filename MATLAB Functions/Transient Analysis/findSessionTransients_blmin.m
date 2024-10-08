function [data] = findSessionTransients_blmin(data,whichstream,whichthreshold,whichblstart,whichblend,whichpostbl,whichcritlabel)
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
% OUTPUTS:
%       DATA:           The original data structure with
%                       sessiontransients_blmin added in.
%                       sessiontransients_blmin includes the max peak
%                       location and value, baseline min location and
%                       value, peak amplitude, and peak rise time in
%                       samples.
%
% Written by R M Donka, June 2024.
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

disp('SESSION TRANSIENTS: Peak baseline determined by min value in the specified baseline window.')

    for eachfile = 1:length(data)
        disp(['Finding Transients: File ',num2str(eachfile)]) % Display which file is being processed
        
        % Find mins and maxes in stream
        allminlocs = find(islocalmin(data(eachfile).(whichstream))); % Find all mins
        allmaxlocs = find(islocalmax(data(eachfile).(whichstream))); % Find all maxes
        if allmaxlocs(1) < allminlocs(1) % If the first max is within a min, remove it
            allmaxlocs = allmaxlocs(2:end);
        end
        
        % Prepare table - preallocate size
        allpkvarnames = {'pkmaxloc', 'pkmaxval', 'pkpreminloc', 'pkpreminval', 'pkamp', 'pkrisesamples', 'pkriseHH_loc', 'pkriseHH_val', 'pkriseHH_samples', 'pkfallHH_loc', 'pkfallHH_val', 'pkfallHH_samples', 'pkHH_AUC'};
        allpkvartypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double','double', 'double', 'double', 'double', 'double', 'double'};
        allpks = table('Size',[length(allmaxlocs), length(allpkvarnames)], 'VariableNames', allpkvarnames, 'VariableTypes', allpkvartypes);

        blstart = data(eachfile).(whichblstart);
        blend = data(eachfile).(whichblend);

        totalnpeaks = 0;
        for eachmax = 1:length(allmaxlocs)
            currpkfallsamples_HH = [];

            currmaxloc = allmaxlocs(eachmax); % Index of current max in data stream
            currmaxval = data(eachfile).(whichstream)(currmaxloc); % Value of current max from data stream

            allpreminlocs = allminlocs(allminlocs > allmaxlocs(eachmax)-blstart & allminlocs < allmaxlocs(eachmax)-blend); % Indexes of all mins that are between the baseline start and end point relative to the current max
            allpreminvals = data(eachfile).(whichstream)(allpreminlocs); % Values of all mins in the window
            currpreminloc = allpreminlocs(find(allpreminvals == min(allpreminvals))); % Index of the baseline window min point in the data stream
            currminval = mean(data(eachfile).(whichstream)(currblstart:currblend)); % Find mean of baseline window

            currpkamp = currmaxval-currminval; % Current peak amplitude

            if currpkamp >= data(eachfile).(whichthreshold) % If the current peak is greater than threshold, add to table allpks
                totalnpeaks = totalnpeaks + 1;
 
                currpktotalrisesamples = currmaxloc - currpreminloc;

                prepkdata = data(eachfile).(whichstream)(currpreminloc:currmaxloc);
                currpkriseval_HH = currminval + currpkamp/2;
                currpkrisesamples_HH = find(prepkdata >= currpkriseval_HH,1,'first');
                currpkriseloc_HH = currpreminloc+currpkrisesamples_HH;

                if (currmaxloc + data(eachfile).(whichpostbl)) > length(data(eachfile).(whichstream))
                    postpkdata = data(eachfile).(whichstream)(currmaxloc:end);
                else
                    postpkdata = data(eachfile).(whichstream)(currmaxloc:(currmaxloc + data(eachfile).(whichpostbl)));
                end
                
                currpkfallsamples_HH = find(postpkdata <= currpkriseval_HH,1,'first');
                currpkfallloc_HH = currmaxloc+currpkfallsamples_HH;
                currpkfallval_HH = data(eachfile).(whichstream)(currpkfallloc_HH);
                              
                currpkdata_HH = data(eachfile).(whichstream)(currpkriseloc_HH:currpkfallloc_HH);
                currpkdata_HH = currpkdata_HH - min(currpkdata_HH);
    
                currauc_HH = round(trapz(currpkdata_HH));

                allpks.pkmaxloc(totalnpeaks) = currmaxloc;
                allpks.pkmaxval(totalnpeaks) = currmaxval;
                allpks.pkpreminloc(totalnpeaks) = currpreminloc;
                allpks.pkpreminval(totalnpeaks) = currminval;
                allpks.pkamp(totalnpeaks) = currpkamp;
                allpks.pkrisesamples(totalnpeaks) = currpktotalrisesamples;
                allpks.pkriseHH_loc(totalnpeaks) = currpkriseloc_HH;
                allpks.pkriseHH_val(totalnpeaks) = currpkriseval_HH;
                allpks.pkriseHH_samples(totalnpeaks) = currpkrisesamples_HH;
                
                if isempty(currpkfallsamples_HH)
                    disp('WARNING: NO FALL FOUND FOR PEAK')
                    disp(totalnpeaks)
                    disp('    Each max index ')
                    disp(eachmax)
                else
                    allpks.pkfallHH_loc(totalnpeaks) = currpkfallloc_HH;
                    allpks.pkfallHH_val(totalnpeaks) = currpkfallval_HH;
                    allpks.pkfallHH_samples(totalnpeaks) = currpkfallsamples_HH;
                    allpks.pkHH_AUC(totalnpeaks) = currauc_HH;
                end
            end               
        end
        data(eachfile).(append('sessiontransients_blmin_',whichcritlabel)) = allpks(1:totalnpeaks,:);
    end
end

