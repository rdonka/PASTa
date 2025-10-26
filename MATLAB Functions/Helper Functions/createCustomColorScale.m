function [customColorScale] = createCustomColorScale(varargin)
% CREATECUSTOMCOLORSCALE  Creates a custom color scale from input hex
% codes and value positions, and outputs color array.
%
%   CREATECUSTOMCOLORSCALE()
%   Creates custom color scale with default paramters.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'hexColors'             - Array; Hex codes of desired colors for color scale.
%                             Default: ["#FFFFFF", "#1CD61C", "#0012FF", "#A600FF", "#FA00FF", "#E00004", "#F26700", "#FFE500"]
%
%   'colorvaluepositions'   - Array; Number positions for colors; Must be
%                             same length as 'hexColors' input.
%                             Default: [0.00, 0.02, 0.183, 0.347, 0.510, 0.673, 0.837, 1.000]
%
%   'noutputcolors'         - Numeric; Number of output colors for the created scale.
%                             Default: 256
%
% OUTPUTS:
%   CUSTOMCOLORSCALE    - Object; Generated color scale for use in heat map
%                         and color grade functions.
%
% EXAMPLE USAGE:
%   customColorScale = createCustomColorScale();
%
%   customColorScale = createCustomColorScale('hexColors',hexcolorarray,'colorvaluepositions',colorvaluepositionarray,...
%                       'noutputcolors', 256);
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Prepare Settings
    % Prepare default values
    defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.

    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'hexColors', defaultparameters.hexColors); % 
    addParameter(p, 'colorvaluepositions', defaultparameters.colorvaluepositions); % 
    addParameter(p, 'noutputcolors', defaultparameters.noutputcolors, @isnumeric); % 

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Display
    disp('CREATECUSTOMCOLORSCALE: Custom color scale created.') % Display file path location
    disp(params)


    hexColors = params.hexColors;
    colorvaluepositions = params.colorvaluepositions;
    noutputcolors = params.noutputcolors;

    
    % Convert HEX -> RGB
    rgbCells  = arrayfun(@(c) sscanf(extractAfter(c,1),'%2x%2x%2x',[1 3]) / 255, hexColors, 'UniformOutput', false);
    rgbColors = vertcat(rgbCells{:});
    
    % Build smooth interpolated map from hexColors and colorvaluepositions
    customColorScale = interp1(colorvaluepositions, rgbColors, linspace(0,1,noutputcolors), 'pchip'); % Linear interpolation
    customColorScale = min(max(customColorScale,0),1); % clamp to [0,1] in case of tiny overshoot