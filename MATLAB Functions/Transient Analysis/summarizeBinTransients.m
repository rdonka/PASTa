function [transientdata] = summarizeBinTransients(transientdata,binfieldname)
% SUMMARIZETRANSIENTS  Assigns each transient to a time bin within the session.
%
%   SUMMARIZETRANSIENTS(TRANSIENTDATA) adds the table
%   'transientsummary_<BINFIELDNAME>' to the TRANSIENTDATA structure, containing
%   frequency and mean quantification of transient events for each session.
%
% REQUIRED INPUTS:
%   TRANSIENTDATA   - Structure array output from findTransients.
%                     Must contain 'params.findTransients' and
%                     'transientquantification' fields.
%
%   BINFIELDNAME    - String; The name of the field in TRANSIENTDATA under
%                     the 'transientquantification' table that contains the
%                     bin IDs for each transient event. For example, 
%                     'Bin_5min'.
%
% OUTPUTS:
%   TRANSIENTDATA   - Original transientdata structure with the summary table
%                     transientsummary_<BINFIELDNAME> added to the specified table of 
%                     transients. Note that the name will change depending
%                     on the BINFIELDNAME specified.
%
% EXAMPLE:
%   transientdata = summarizeTransients(transientdata);
%
% See also: findTransients
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Display
    disp('SUMMARIZEBINTRANSIENTS: Summarizing session transients by bin for each session. Summary data will be added to TRANSIENTDATA as transientsummary_session.')
    disp('   Output table includes each bin as a row with transient frequency and session means of all other transient quantification variables.')
    disp(['     binfieldname: ',binfieldname])

    %% Summarize Transients
    for eachfile = 1:length(transientdata)
        % Prepare variables
        fs = transientdata(eachfile).params.findTransients.fs;
        binsamples = transientdata(eachfile).params.binTransients.(binfieldname).binlengthsamples;
        binlengthmins = transientdata(eachfile).params.binTransients.(binfieldname).binlengthmins;
        binlengthseconds = binlengthmins*60;

        uniquebins = unique(transientdata(eachfile).transientquantification.(binfieldname));

        % Display number of bins:
        disp(['File ',num2str(eachfile),': ',num2str(length(uniquebins)),' bins identified.'])
        % Summarize transient means by subject for whole session
        transientsummary_bin = table(uniquebins, 'VariableNames',{binfieldname});
        
        for eachbin = 1:length(uniquebins)
            currbinval = uniquebins(eachbin);
            binlocs = find(transientdata(eachfile).transientquantification.(binfieldname) == currbinval);

            transientsummary_bin.freq(eachbin) = height(transientdata(eachfile).transientquantification.transientID(binlocs));
            transientsummary_bin.freqpermin(eachbin) = height(transientdata(eachfile).transientquantification.transientID(binlocs))/binlengthmins;
            transientsummary_bin.freqhz(eachbin) = height(transientdata(eachfile).transientquantification.transientID(binlocs))/binlengthseconds;
    
            transientsummary_bin.maxval(eachbin) = mean(transientdata(eachfile).transientquantification.maxval(binlocs), 'omitmissing');
            transientsummary_bin.blval(eachbin) = mean(transientdata(eachfile).transientquantification.blval(binlocs), 'omitmissing');
            transientsummary_bin.amp(eachbin) = mean(transientdata(eachfile).transientquantification.amp(binlocs), 'omitmissing');
            transientsummary_bin.quantheightval(eachbin) = mean(transientdata(eachfile).transientquantification.quantheightval(binlocs), 'omitmissing');
            transientsummary_bin.risesamples(eachbin) = mean(transientdata(eachfile).transientquantification.risesamples(binlocs), 'omitmissing');
            transientsummary_bin.risems(eachbin) = mean(transientdata(eachfile).transientquantification.risems(binlocs), 'omitmissing');
            transientsummary_bin.fallsamples(eachbin) = mean(transientdata(eachfile).transientquantification.fallsamples(binlocs), 'omitmissing');
            transientsummary_bin.fallms(eachbin) = mean(transientdata(eachfile).transientquantification.fallms(binlocs), 'omitmissing');
            transientsummary_bin.widthsamples(eachbin) = mean(transientdata(eachfile).transientquantification.widthsamples(binlocs), 'omitmissing');
            transientsummary_bin.widthms(eachbin) = mean(transientdata(eachfile).transientquantification.widthms(binlocs), 'omitmissing');
            transientsummary_bin.AUC(eachbin) = mean(transientdata(eachfile).transientquantification.AUC(binlocs), 'omitmissing');
            transientsummary_bin.IEIsamples(eachbin) = mean(transientdata(eachfile).transientquantification.IEIsamples(binlocs), 'omitmissing');        
            transientsummary_bin.IEIms(eachbin) = mean(transientdata(eachfile).transientquantification.IEIms(binlocs), 'omitmissing');
            transientsummary_bin.IEIs(eachbin) = mean(transientdata(eachfile).transientquantification.IEIs(binlocs), 'omitmissing');
            transientsummary_bin.compoundeventtotal(eachbin) = sum(transientdata(eachfile).transientquantification.compoundeventnum(binlocs) > 0);
        end
        transientdata(eachfile).(append('transientsummary_',binfieldname)) = transientsummary_bin;
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