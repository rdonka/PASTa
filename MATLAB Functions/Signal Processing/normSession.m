function [data] = normSession(data,whichstream)
% NORMSESSION    Normalizes the data stream to zscore based on the whole session.
% INPUTS:
%       DATA:           This is a data structure that contains at least the
%                       stream specified to be normalized.
%
%       WHICHSTREAM:    The name of the field containing the stream to be
%                       normalized.
% OUTPUTS:
%       DATA:           The original data structure with data.sigz_normsession added.

% Written by R M Donka, August 2024.
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

    for eachfile = 1:length(data)
        data(eachfile).sigz_normsession = zscore(data(eachfile).(whichstream));
    end
end
