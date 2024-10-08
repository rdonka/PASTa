function [data] = binSessionTransients(data,whichbinsamples,whichstream,whichtransients,whichpklocs)
% FINDSESSIONTRANSIENTS   Finds transients for the whole session of data.
%
% INPUTS:
%       DATA:               This is a structure that contains at least a
%                           field containing the transients you want to bin.
%
%       WHICHBINSAMPLES:    The name of the field containing the number of
%                           samples per bin.
%
%       WHICHSTREAM:        The name of the field containing the data
%                           stream input to the findSessionTransients to
%                           identify transients from.
%
%       WHICHTRANSIENTS:    The name of the field containing the transients 
%                           that you want to identify bins for.
%
%       WHICHPKLOCS:        The name of the field containing the transient
%                           peak locations (indexes) relative to the whole
%                           session.
%
% OUTPUTS:
%       DATA:               This is the original data structure with bins
%                           added to the specified table of transients
%
% Written by R M Donka, March 2024.
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.
    for eachfile = 1:length(data)
        disp(['Binning Transients: File ',num2str(eachfile)]) % Display which file is being processed
        
        binsamples = data(eachfile).(whichbinsamples);

        % Bin the transients into 5 minute bins
        nbins = ceil(length(data(eachfile).(whichstream))/binsamples);
        disp(['     Total number of bins: ',num2str(nbins)]) % Display which file is being processed

        for eachbin = 0:(nbins-1)
            startbin = (eachbin*binsamples)+1;
            endbin = (eachbin+1)*binsamples;
            data(eachfile).(whichtransients).Bin((data(eachfile).(whichtransients).(whichpklocs) >= startbin & ...
               data(eachfile).(whichtransients).(whichpklocs) < endbin)) = eachbin+1;
        end
    end

end