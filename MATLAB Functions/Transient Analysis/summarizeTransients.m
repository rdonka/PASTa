function [transientdata] = summarizeTransients(transientdata)
% SUMMARIZETRANSIENTS  Assigns each transient to a time bin within the session.
%
%   SUMMARIZETRANSIENTS(TRANSIENTDATA) adds the table
%   'transientsummary_session' to the TRANSIENTDATA structure, containing
%   frequency and mean quantification of transient events for each session.
%
% REQUIRED INPUTS:
%   TRANSIENTDATA   - Structure array output from findTransients.
%                     Must contain 'params.findTransients' and 
%                     'transientquantification' fields.
%
% OUTPUTS:
%   TRANSIENTDATA   - Original transientdata structure with the summary table
%                     transientsummary_session added to the specified table of 
%                     transients.
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
    disp('SUMMARIZETRANSIENTS: Summarizing session transients for each session. Summary data will be added to TRANSIENTDATA as transientsummary_session.')
    disp('   Output table includes transient frequency and session means of all other transient quantification variables.')

    %% Summarize Transients
    for eachfile = 1:length(transientdata)
        % Prepare variables
        fs = transientdata(eachfile).params.findTransients.fs;
        streamtotalsamples = transientdata(eachfile).params.findTransients.streamtotalsamples;
        streamtotalseconds = streamtotalsamples/fs;
        streamtotalminutes = streamtotalsamples/fs/60;

        % Summarize transient means by subject for whole session
        transientsummary_session = table();
        transientsummary_session.freq = height(transientdata(eachfile).transientquantification);
        transientsummary_session.freqpermin = height(transientdata(eachfile).transientquantification)/streamtotalminutes;
        transientsummary_session.freqhz = height(transientdata(eachfile).transientquantification)/streamtotalseconds;

        transientsummary_session.maxval = mean(transientdata(eachfile).transientquantification.maxval, 'omitmissing');
        transientsummary_session.blval = mean(transientdata(eachfile).transientquantification.blval, 'omitmissing');
        transientsummary_session.amp = mean(transientdata(eachfile).transientquantification.amp, 'omitmissing');
        transientsummary_session.quantheightval = mean(transientdata(eachfile).transientquantification.quantheightval, 'omitmissing');
        transientsummary_session.risesamples = mean(transientdata(eachfile).transientquantification.risesamples, 'omitmissing');
        transientsummary_session.risems = mean(transientdata(eachfile).transientquantification.risems, 'omitmissing');
        transientsummary_session.fallsamples = mean(transientdata(eachfile).transientquantification.fallsamples, 'omitmissing');
        transientsummary_session.fallms = mean(transientdata(eachfile).transientquantification.fallms, 'omitmissing');
        transientsummary_session.widthsamples = mean(transientdata(eachfile).transientquantification.widthsamples, 'omitmissing');
        transientsummary_session.widthms = mean(transientdata(eachfile).transientquantification.widthms, 'omitmissing');
        transientsummary_session.AUC = mean(transientdata(eachfile).transientquantification.AUC, 'omitmissing');
        transientsummary_session.IEIsamples = mean(transientdata(eachfile).transientquantification.IEIsamples, 'omitmissing');        
        transientsummary_session.IEIms = mean(transientdata(eachfile).transientquantification.IEIms, 'omitmissing');
        transientsummary_session.IEIs = mean(transientdata(eachfile).transientquantification.IEIs, 'omitmissing');
        transientsummary_session.compoundeventtotal = sum(transientdata(eachfile).transientquantification.compoundeventnum > 0);

        transientdata(eachfile).transientsummary_session = transientsummary_session;
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