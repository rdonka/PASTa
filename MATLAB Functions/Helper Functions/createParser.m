function p = createParser(functionName)
% CREATEPARSER    Returns an inputParser object with universal settings.
%                 Requires input names to be fully matched and excludes
%                 any that do not match specified inputs.
%
% Copyright (C) 2025 Rachel Donka. Licensed under the GNU General Public License v3.

    % Create a blank parser
    p = inputParser();

    % Set the name of the function the parser is for
    if nargin < 1 || isempty(functionName)
        p.FunctionName = 'MyFunction';
    else
        p.FunctionName = functionName;
    end

    % Set universal defaults
    p.KeepUnmatched    = false;
    p.PartialMatching  = false;
end