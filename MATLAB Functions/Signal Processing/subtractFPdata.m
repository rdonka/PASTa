function [data] = subtractFPdata(data, whichsigfield, whichbaqfield, whichfs, varargin)
% SUBTRACTFPDATA    Subtracts background fluorescence from signal and applies
%                   filtering to fiber photometry data.
%
%   DATA = SUBTRACTFPDATA(DATA, WHICHSIGFIELD, WHICHBAQFIELD, WHICHFS, 'PARAM1', VAL1, ...)
%   processes the fluorescence data by subtracting the specified background
%   channel from the signal channel and applying a Butterworth filter. The
%   function outputs the data as either delta F/F ('dff') or delta F ('df').
%
%
% REQUIRED INPUTS:
%       DATA            - Struct array; each element contains fields for
%                         signal and background fluorescence data, as well
%                         as the sampling rate.
%
%       WHICHSIGFIELD   - String; name of the field in DATA containing the
%                         "signal" data stream (e.g., 'sig').
%
%       WHICHBAQFIELD   - String; name of the field in DATA containing the
%                         "background" data stream (e.g., 'baq').
%
%       WHICHFS         - String; name of the field in DATA containing the
%                         sampling rate (e.g., 'fs').
%
% PTIONAL INPUT NAME-VALUE PAIRS:
%       'baqscalingtype'     - String; method for scaling the background
%                              fluorescence before subtraction. Options are:
%                              'frequency' (default), 'sigmean', 'OLS',
%                              'detrendedOLS', 'smoothedOLS', 'IRLS'. See
%                              PASTa user guide signal processing section
%                              for more details.
%
%       'baqscalingfreqmin'  - Numeric; minimum frequency (Hz) for scaling
%                              when 'frequency' method is used. Default: 10
%
%       'baqscalingfreqmax'  - Numeric; maximum frequency (Hz) for scaling
%                              when 'frequency' method is used. Default: 100
%
%       'baqscalingperc'     - Numeric; percentage to adjust the background
%                              scaling factor. Default: 1 (100%)
%
%       'subtractionoutput'  - String; specifies the output type:
%                              'dff' (default) for delta F/F, 'df' for delta F
%
%       'artifactremoval'    - Logical (0 or 1); set to 1 to detect and remove
%                              artifacts from data streams. Default: 0 (false)
%
%       'filtertype'         - String; Type of filter to apply:
%                              'bandpass' (default), 'highpass', 'lowpass',
%                              'nofilter'.
%
%       'padding'            - Logical (0 or 1); set to 1 to apply padding 
%                              before filtering. Default: 1 (true)
%
%       'paddingperc'        - Numeric; percentage of data length to use for
%                              padding. Minimum 0.1 (10%). Default: 0.1
%
%       'filterorder'        - Numeric; order of the Butterworth filter.
%                              Default: 3
%
%       'highpasscutoff'     - Numeric; high-pass filter cutoff frequency (Hz).
%                              Default: 0.0051
%
%       'lowpasscutoff'      - Numeric; low-pass filter cutoff frequency (Hz).
%                              Default: 2.2860
%
%       'suppressdisp'       - Logical (0 or 1); set to 1 to suppress 
%                              command window displays. Default: 0 (false)
%
%   OUTPUT:
%       DATA   - Struct array; the input DATA structure with additional fields:
%              - 'baqscaled': scaled background fluorescence data.
%              - 'sigsub': subtracted signal (delta F or delta F/F).
%              - 'sigfilt': filtered subtracted signal.
%              - 'inputs_subtractFPdata': structure containing input parameters.
%              - Additional fields related to artifact removal if enabled.
%
%   EXAMPLES:
%       % Example 1: Basic subtraction with default settings
%       data = subtractFPdata(data, 'sig', 'baq', 'fs');
%
%       % Example 2: Subtraction with specified background scaling and filter
%       data = subtractFPdata(data, 'sig', 'baq', 'fs', ...
%                             'baqscalingtype', 'IRLS', ...
%                             'filtertype', 'lowpass', ...
%                             'lowpasscutoff', 1.5);
%
%   NOTES:
%       - Ensure that the signal and background fields in DATA have the same
%         length. If they differ, the function will adjust accordingly and
%         issue a warning.
%       - The function adds new fields to the DATA structure to store results
%         and input parameters.
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

    % Add optional name-value pair arguments with validation
    addParameter(p, 'baqscalingtype', defaultparameters.baqscalingtype, @(x) ischar(x) && ismember(x, {'frequency', 'sigmean', 'OLS', 'detrendedOLS', 'smoothedOLS', 'IRLS'})); % baqscalingtype: input subtraction type must be one of 6 options
    addParameter(p, 'baqscalingfreqmin', defaultparameters.baqscalingfreqmin, @(x) validateattributes(x, {'numeric'}, {'positive'})); % baqscalingfreqmin: input must be numeric and positive
    addParameter(p, 'baqscalingfreqmax', defaultparameters.baqscalingfreqmax, @(x) validateattributes(x, {'numeric'}, {'positive'})); % baqscalingfreqmax: input must be numeric and positive
    addParameter(p, 'baqscalingperc', defaultparameters.baqscalingperc, @(x) validateattributes(x, {'numeric'}, {'positive'})); % baqscalingperc: input must be numeric and positive
    addParameter(p, 'subtractionoutput', defaultparameters.subtractionoutput, @(x) ischar(x) && ismember(x, {'dff', 'df'})); % subtractionoutput: input must be either 'dff' or 'df'
    addParameter(p, 'artifactremoval', defaultparameters.artifactremoval, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1]))); % artifactremoval: input must be logical or numeric (either 0 or 1)
    addParameter(p, 'filtertype', defaultparameters.filtertype, @(x) ischar(x) && ismember(x, {'nofilter', 'bandpass', 'highpass', 'lowpass'})); % filtertype: input must be one of 4 options
    addParameter(p, 'padding', defaultparameters.padding, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1]))); % padding: input must be logical or numeric (either 0 or 1)
    addParameter(p, 'paddingperc', defaultparameters.paddingperc, @(x) validateattributes(x, {'numeric'}, {'positive', '<=', 1,'>=',0.1})); % paddingperc: input must be between 10% and 100%
    addParameter(p, 'filterorder', defaultparameters.filterorder, @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'})); % filterorder: input must be a positive integer
    addParameter(p, 'highpasscutoff', defaultparameters.highpasscutoff, @(x) validateattributes(x, {'numeric'}, {'positive'})); % highpasscutoff: input must be numeric and positive
    addParameter(p, 'lowpasscutoff', defaultparameters.lowpasscutoff, @(x) validateattributes(x, {'numeric'}, {'positive'})); % lowpasscutoff: input must be numeric and positive
    addParameter(p, 'suppressdisp', defaultparameters.suppressdisp, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1]))); % suppressdisp: input must be logical or numeric (either 0 or 1)

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Main display and function inputs
    if params.suppressdisp == 0
        if strcmp(params.baqscalingtype,'frequency')==true % Display for frequency background scaling - default
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with frequency domain scaling (threshold at ', num2str(params.baqscalingfreqmin), 'hz). Subtracted signal will be output as ', params.subtractionoutput, ' to data.sigsub.'));
        elseif strcmp(params.baqscalingtype,'sigmean')==true % Display for signal mean background scaling
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with signal mean scaling. Subtracted signal will be output as ', params.subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(params.baqscalingtype,'OLS')==true % Display for ordinary least squares regression background scaling
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with ordinary least-squares regression (OLS). Subtracted signal will be output as ', params.subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(params.baqscalingtype,'detrendedOLS')==true % Display for linear detrending and ordinary least squares regression background scaling
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with detrending and ordinary least-squares regression (OLS). Subtracted signal will be output as ', params.subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(params.baqscalingtype,'smoothedOLS')==true % Display for time domain subtraction
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with lowess smoothing and ordinary least-squares regression (OLS). Subtracted signal will be output as ', params.subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(params.baqscalingtype,'IRLS')==true % Display for time domain subtraction
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with iteratively reweighted least squares regression (IRLS). Subtracted signal will be output as ', params.subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        else
            disp(append('ERROR: Baq scaling type issue - baqscalingtype set to ', params.baqscalingtype, '. Function inputs limited to frequency, sigmean, OLD, detrendedOLS, smoothedOLS, or IRLS.'));
        end

        if strcmp(params.filtertype,'nofilter')==true % Display for filter settings
            disp(append('WARNING: filtertype manually set to nofilter - NO FILTER APPLIED'));
        else
            disp(append('Filter type set to ', params.filtertype, '. Subtracted and filtered signal will be output as data.sigfilt.'));
            if params.padding == 1 % Padding application display
                disp('   NOTE: Padding applied prior to filtering and removed from final output.')
            end
        end
        if params.artifactremoval == 1 % Display if artifact removal is turned on
            disp(append('Artifact removal set to 1 - artifacts will be automatically detected and replaced with NaNs in subtracted and filtered outputs.'));
        end

        disp('   PARAMETERS:') % Display all input values
        disp(params)

        disp('    Parameters using default values:'); % Display input parameters values set to default values
        disp(p.UsingDefaults);
    end

%% Subtract and Filter Data
    for eachfile = 1:length(data)
        data(eachfile).inputs_subtractFPdata = params; % Add inputs to data frame

        if params.suppressdisp == 0
            fprintf('Subtracting file number: %.f \n',eachfile) % Display which file is being subtracted
        end
        
        % Prepare sampling rate
        fs = data(eachfile).(whichfs);
        
        % Prepare filters
        if strcmp(params.filtertype,'bandpass')==true
            highpasscutoffval = round(params.highpasscutoff/floor(fs/2),6);
            lowpasscutoffval = round(params.lowpasscutoff/floor(fs/2),6);
            [a,b] = butter(params.filterorder,highpasscutoffval,'high');
            [c,d] = butter(params.filterorder,lowpasscutoffval,'low');
        elseif strcmp(params.filtertype, 'highpass')==true
            highpasscutoffval = round(params.highpasscutoff/floor(fs/2),6);
            [a,b] = butter(params.filterorder,highpasscutoffval,'high');
        elseif strcmp(params.filtertype, 'lowpass')==true
            lowpasscutoffval = round(params.lowpasscutoff/floor(fs/2),6);
            [c,d] = butter(params.filterorder,lowpasscutoffval,'low');
        end

        try
        %% Prepare data: Check length match of signal and background
            % Check length of streams and make temp stream variables
            if length(data(eachfile).(whichsigfield))==length(data(eachfile).(whichbaqfield)) % If the length of 405 and 465 are equal, set temps to raw values
                baq = data(eachfile).(whichbaqfield);
                sig = data(eachfile).(whichsigfield);
            elseif length(data(eachfile).(whichsigfield))>length(data(eachfile).(whichbaqfield)) % If the length of signal is greater, set temps to length of 405
                disp('   WARNING: Length of signal greater than length of background stream. Signal data adjusted to match background length.')
                baq = data(eachfile).(whichbaqfield);
                sig = data(eachfile).(whichsigfield)(1:length(data(eachfile).(whichbaqfield)));
            elseif length(data(eachfile).(whichsigfield))<length(data(eachfile).(whichbaqfield)) % If the length of 405 is greater, set temps to length of 465
                disp('   WARNING: Length of background greater than length of signal stream. Background data adjusted to match signal length.')
                sig = data(eachfile).(whichsigfield);
                baq = data(eachfile).(whichbaqfield)(1:length(data(eachfile).(whichsigfield)));
            end

        %% Scale Background: apply the selected background scaling method.
            if strcmp('frequency',params.baqscalingtype)==true % Frequency domain scaling - scales to power
                baq_centered = baq - mean(data(eachfile).(whichbaqfield)); % Center background at 0
                sig_centered = sig - mean(data(eachfile).(whichsigfield)); % Center signal at 0

                [sigFFT, sigF] = preparestreamFFT(sig_centered,fs); % Prep FFT
                [baqFFT, baqF] = preparestreamFFT(baq_centered,fs); % Prep FFT

                sigFidxs = sigF > params.baqscalingfreqmin & sigF < params.baqscalingfreqmax; % Find indices of frequencies above the set min threshold and below the set max threshold
                baqFidxs = baqF > params.baqscalingfreqmin & baqF < params.baqscalingfreqmax; % Find indices of frequencies above the set min threshold and below the set max threshold

                sigPower = sigFFT(sigFidxs).^2; % Compute average power in the selected band
                baqPower = baqFFT(baqFidxs).^2; % Compute average power in the selected band

                powerRatio = mean(sigPower) / mean(baqPower); % Find the ratio of mean powers between signal and background (same as RMS)
                baqscalingfactor = sqrt(powerRatio) * params.baqscalingperc;  % or remove baqscalingperc if not needed

                baqscaled = (baq_centered*baqscalingfactor) + mean(sig); % Adjust back to same units as raw signal
                data(eachfile).baqscalingfactor = baqscalingfactor; % Add constant baqscalingfactor to data

                if baqscalingfactor > 3 % Display a warning if the baqscaling factor is high
                    disp(append('   WARNING: Possible over scaling. baqscalingfactor = ',num2str(baqscalingfactor)))
                end
            elseif strcmp('frequency_amplitude',params.baqscalingtype)==true % Frequency domain scaling - scales to amplitude
                baq_centered = baq - mean(data(eachfile).(whichbaqfield)); % Center background at 0
                sig_centered = sig - mean(data(eachfile).(whichsigfield)); % Center signal at 0

                [sigFFT, sigF] = preparestreamFFT(sig_centered,fs); % Prep FFT
                [baqFFT, baqF] = preparestreamFFT(baq_centered,fs); % Prep FFT

                sigFidxs = sigF > params.baqscalingfreqmin; % Find indices of frequencies above the set threshold
                baqFidxs = baqF > params.baqscalingfreqmin; % Find indices of frequencies above the set threshold

                baqscalingfactor = (mean(sigFFT(sigFidxs))/mean(baqFFT(baqFidxs)))*params.baqscalingperc; % Find power ratio of signal to background for frequencies above the set threshold
                baqscaled = (baq_centered*baqscalingfactor) + mean(sig); % Adjust back to same units as raw signal

                data(eachfile).baqscalingfactor = baqscalingfactor; % Add constant baqscalingfactor to data
            elseif strcmp('sigmean',params.baqscalingtype)==true % Signal mean time domain scaling
                baqscalingfactor = (mean(sig)/(mean(baq)))*params.baqscalingperc; % Find the scaling factor for background to signal based on the ratio of the means
                baqscaled = baq*baqscalingfactor;
                data(eachfile).baqscalingfactor = baqscalingfactor; % Add constant baqscalingfactor to data
            elseif strcmp('OLS',params.baqscalingtype)==true % OLS regression scaling; Per GuPPY, see Sherathiya et al 2021, https://www.nature.com/articles/s41598-021-03626-9
                bls=polyfit(baq,sig,1);
                baqscaled=bls(1).*baq+bls(2);      
            elseif strcmp('detrendedOLS',params.baqscalingtype)==true % Linear detrend and OLS regression scaling
                sigmean = mean(sig);
                baq = detrend(baq); % Detrend background
                sig = detrend(sig)+sigmean; % Detrend signal
                bls = polyfit(baq,sig,1); % Fit background to signal with least-squares linear regression
                baqscaled = (bls(1).*baq)+ bls(2); % Fit to signal
                data(eachfile).baqdetrend = baq; % Add detrended signal to data
                data(eachfile).sigdetrend = sig; % Add detrended background to data
            elseif strcmp('smoothedOLS',params.baqscalingtype)==true % Lowess smoothing and OLS regression scaling; Per pMAT, see Bruno et al 2021, https://www.sciencedirect.com/science/article/abs/pii/S0091305720307413?via=ihub#f0020
                baq=smooth(baq,0.002,'lowess')'; 
                sig=smooth(sig,0.002,'lowess')';
                bls=polyfit(baq(1:end),sig(1:end),1);
                baqscaled=bls(1).*baq+bls(2);
                data(eachfile).baqsmoothed = baq; % Add detrended signal to data
                data(eachfile).sigsmoothed = sig; % Add detrended background to data
            elseif strcmp('IRLS',params.baqscalingtype)==true % IRLS regression scaling; See Keevers et al 2024, https://www.researchsquare.com/article/rs-3549461/v2
                IRLS_coeffs = reshape(flipud(robustfit(baq, sig, 'bisquare', 1.4, 'on')), [1, 2]);
                baqscaled = polyval(IRLS_coeffs,baq);
            else
                disp(append('ERROR: Baq scaling type issue - baqscalingtype set to ', params.baqscalingtype, '. Function inputs limited to frequency, sigmean, OLS, detrendedOLS, smoothedOLS, or IRLS.'));
            end
            data(eachfile).(append(whichbaqfield,'scaled')) =  baqscaled; % Add scaled background to data frame

        %% Subtract data
            if strcmp(params.subtractionoutput,'dff')==true % Output delta F/F
                data(eachfile).sigsub = (sig-baqscaled)./baqscaled*100;
            elseif strcmp(params.subtractionputput,'df')==true % Output delta F
                data(eachfile).sigsub = (sig-baqscaled);
            else
                disp(append('ERROR: Subtraction Output type issue - subtractionoutput set to ', subtractionoutput, '. Function inputs limited to df or dff.'));
            end
           
        %% OPTIONAL: Remove Artifacts
            if params.artifactremoval == 1
                artifactremovaldata = removeStreamArtifacts([data(eachfile).sigsub],'sigsub',fs);
                if artifactremovaldata.numartifacts > 0
                    disp(append('   WARNING File ',num2str(eachfile),': ',num2str(artifactremovaldata.numartifacts),' artifacts found.'))
                end
                data(eachfile).sigsubraw = data(eachfile).sigsub; % Raw unadjusted sigsub is added as sigsub_raw
                data(eachfile).sigsubAR = artifactremovaldata.sigsubAR; % sigsub with artifacts replaced the NaNs is added as a separate field
                data(eachfile).sigsub = artifactremovaldata.sigsubAM; % sigsub is replaced with the sigsub with artifacts replaced with signal mean
                data(eachfile).numartifacts = artifactremovaldata.numartifacts; % Add number of artifacts   
                data(eachfile).artifacts = artifactremovaldata.artifacts; % Add artifact indexes and values
            end

        %% Filter subtracted signal
            if strcmp(params.filtertype, 'nofilter')==true % Skip filtering
            else
                if params.padding==1 % Filter with padding
                    nsamplesedge = floor(length(data(eachfile).sigsub)*params.paddingperc); % Padding: determine number of samples to append to start and end
                    firstsamples = [flip(data(eachfile).sigsub(1:nsamplesedge),2)]; % Padding: Extract data to append to beginning of signal
                    lastsamples = [flip(data(eachfile).sigsub(length(data(eachfile).sigsub)-(nsamplesedge-1):length(data(eachfile).sigsub)),2)]; % Padding: Extract data to append to end of signal
                    sigsubdata = [firstsamples data(eachfile).sigsub lastsamples]; % Padding: Create padded data
                    
                    % Apply filter
                    if strcmp(params.filtertype,'bandpass') 
                        filtdata = filtfilt(a,b,sigsubdata); % Apply the highpass
                        filtdata = filtfilt(c,d,filtdata); % Apply the lowpass
                    elseif strcmp(params.filtertype,'highpass')
                        filtdata = filtfilt(a,b,sigsubdata); % Apply the highpass
                    elseif strcmp(params.filtertype,'lowpass')
                        filtdata = filtfilt(c,d,sigsubdata); % Apply the highpass
                    else
                        disp('WARNING: Filter type undefined - no filter applied.')
                    end
                    data(eachfile).sigfilt = filtdata(nsamplesedge+1:length(filtdata)-nsamplesedge); % Padding: Remove padding from beginning and end
                else % Filter without padding
                    if strcmp(params.filtertype,'bandpass')
                        filtdata = filtfilt(a,b,data(eachfile).sigsub); % Apply the highpass
                        filtdata = filtfilt(c,d,filtdata); % Apply the lowpass
                    elseif strcmp(params.filtertype,'highpass')
                        filtdata = filtfilt(a,b,data(eachfile).sigsub); % Apply the highpass
                    elseif strcmp(params.filtertype,'lowpass')
                        filtdata = filtfilt(c,d,data(eachfile).sigsub); % Apply the highpass
                    else
                        disp('WARNING: Filter type undefined - no filter applied.')
                    end
                    data(eachfile).sigfilt = filtdata;
                end
            end

            %% OPTIONAL: Remove artifact periods from filtered signal
            % Filtering can't include NaNs in streams, so sigsub with artifacts replaced by signal mean is passed to filter. 
            % Filtered stream is adjusted here to replace artifact periods with NaNs for future analyses.
            if params.artifactremoval == 1
                nansigfilt = data(eachfile).sigfilt; % Temp array of filtered signal with artifacts replaced by signal mean
                if artifactremovaldata.numartifacts > 0
                    for eachartifact = 1:data(eachfile).numartifacts % Remove each artifact from stream and replace with NaNs
                        nansigfilt(data(eachfile).artifacts.startloc(eachartifact):data(eachfile).artifacts.endloc(eachartifact)) = NaN; 
                    end
                end
                data(eachfile).sigfiltAR = nansigfilt; % Add fitlered signal with artifacts replaced with NaNs to the data structure
            end

        catch
        disp(['WARNING: subtractFPdata without success - file ' num2str(eachfile)]); % Error message to display
        continue;
        end
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