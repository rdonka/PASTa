function [data] = normBaseline(data,whichstream,whichblstart,whichblend)
% NORMSESSION    Normalizes the whole session data stream based on a
%                baseline period.
% INPUTS:
%       DATA:           Data structure; Must contain at least the stream to 
%                       be normalized, baseline start field, and baseline 
%                       end field.
%
%       WHICHSTREAM:    String; The name of the field containing the data 
%                       stream to be normalized.
%
%       WHICHBLSTART:   String; The name of the field containing the index 
%                       of the start of the baseline period.
%
%       WHICHBLEND:     String; The name of the field containing the index 
%                       of the end of the baseline period.
%
% OUTPUTS:
%       DATA:           Data structure; This is the original data structure 
%                       with 'data.WHICHSTREAMz_normbaseline' added.
%
% Written by R M Donka, August 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Normalize to session baseline
disp(append('NORM BASELINE: Normalizing ',whichstream,' to session baseline mean and standard deviation.'))
disp(append('     Baseline defined by ',whichblstart,' and ',whichblend,'.'))

    for eachfile = 1:length(data)
        disp(append('   Normalizing: File ',num2str(eachfile)))
        try
            BLmean = mean(data(eachfile).(whichstream)(data(eachfile).(whichblstart):data(eachfile).(whichblend))); % Find the mean of the session baseline
            BLsd = std(data(eachfile).(whichstream)(data(eachfile).(whichblstart):data(eachfile).(whichblend))); % Find the standard deviation of the session baseline
    
            data(eachfile).(append(whichstream,'z_normbaseline')) = (data(eachfile).(whichstream) - BLmean)/BLsd; % Z score the whole session to the baseline
        catch
            disp(append('WARNING: File ',num2str(eachfile), ' - failed to normalize stream: ', whichstream)) 
        end
    end
end
