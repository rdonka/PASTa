function [data] = normSession(data,whichstream)
% NORMSESSION    Normalizes the data stream to zscore based on the whole session.
% INPUTS:
%       DATA:           Data structure; Must contain at least the stream 
%                       specified to be normalized.
%
%       WHICHSTREAM:    String; The name of the field containing the stream 
%                       to be normalized.
% OUTPUTS:
%       DATA:           Data structure; The original data structure with 
%                       'data.WHICHSTREAMz_normsession' added.
%
% Written by R M Donka, August 2024.
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Normalize to whole session
disp(append('NORM SESSION: Normalizing ',whichstream,' to whole session mean and standard deviation.'))

    for eachfile = 1:length(data)
        disp(append('   Normalizing: File ',num2str(eachfile)))
        data(eachfile).(append(whichstream, 'z_normsession')) = zscore(data(eachfile).(whichstream)); % Z score whole session
    end
end
