function [data] = normBaseline(data,whichstream,BLstart,BLend)
% NORMSESSION    Normalizes the whole session data stream based on a
%                baseline period.
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       stream to be normalized, baseline start, and 
%                       baseline end sample numbers.
%
%       NORMTYPE:       This is a character string to specify the type of
%                       normalization to apply.
%                       'df': delta f/f - Normalizes the session to the max
%                             value of the signal.
%                       'z': z score - Normalizes the raw values based on
%                            the whole session mean and standard deviation.
%                       'both': Returns both df/f and z scored data.  
%
%       WHICHSTREAM:    The name of the field containing the data stream to
%                       be normalized.
%
%       BLSTART:        The name of the field in the data structure to use
%                       for the start of the baseline period.
%
%       BLEND:          The name of the field in the data structure to use
%                       for the end of the baseline period.
%
% OUTPUTS:
%       DATA:           This is the original data structure with 
%                       data.(streamname)_normbaseline_z added.

% Written by R M Donka, January 2024.
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes
    for eachfile = 1:length(data)
        disp(append('   NORMALIZING FILE NUMBER ',num2str(eachfile)))

        try
            BLmean = mean(data(eachfile).(whichstream)(data(eachfile).(BLstart):data(eachfile).(BLend)));
            BLsd = std(data(eachfile).(whichstream)(data(eachfile).(BLstart):data(eachfile).(BLend)));
    
            data(eachfile).(append(whichstream,'sigz_normbaseline')) = (data(eachfile).(whichstream) - BLmean)/BLsd;
        catch
            disp(append('WARNING: File number ',num2str(eachfile), ' - failed to normalize stream: ', whichstream))
        end
    end
end
