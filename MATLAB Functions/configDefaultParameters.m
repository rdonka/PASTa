function defaultparameters = configDefaultParameters(callerFunction)
% CONFIGDEFAULTPARAMETERS  Returns a struct of default parameters for various functions.
%
%   DEFAULTPARAMETERS = CONFIGDEFAULTPARAMETERS(CALLERFUNCTION) returns only
%   the parameters relevant to the specified function. If no function name is
%   provided, all parameters are returned.
%
% INPUTS:
%   CALLERFUNCTION - (optional) String; Name of the function requesting parameters.
%                    Example: 'cropFPdata', 'extractTDTdata'
%
% OUTPUTS:
%   DEFAULTPARAMETERS - Struct containing relevant parameters.
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

% Define all default parameters for different functions
    % extractTDTdata
    allparameters.extractTDTdata.trim = 5;
    allparameters.extractTDTdata.skipexisting = 1;

    % subtractFPdata




    % removeStreamArtifacts
    allparameters.removeStreamArtifacts.outlierthresholdk = 10; % Outlier detection threshold; k = 3 or 4 is a common choice for outlier detection
    allparameters.removeStreamArtifacts.artifactremovalwindow = 0.3; % Seconds to replace with NaNs before and after artifact
    allparameters.removeStreamArtifacts.artifactampthreshold_max = 8; % Artifact SD detection threshold - max artifacts
    allparameters.removeStreamArtifacts.artifactampthreshold_min = 8; % Artifact SD detection threshold - min artifacts
    allparameters.removeStreamArtifacts.bucketsizeSecs = 30; % Bucket size (seconds) for session amplitude and mean calculations



% If a specific function is requested, return only its parameters
    if nargin > 0 && isfield(allparameters, callerFunction)
        defaultparameters = allparameters.(callerFunction);
    else
        % If no function is specified, return the entire structure
        defaultparameters = allparameters;
    end
end
