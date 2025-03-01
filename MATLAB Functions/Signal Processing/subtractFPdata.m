function [data] = subtractFPdata(data,whichsigfield,whichbaqfield,whichfs,varargin)
% SUBTRACTFPDATA    Subtracts the specified background channel from the 
%                   specified signal channel and applies a Butterworth 
%                   filter. Outputs data as dF/F (default) or dF.
%
% Copyright (C) 2025 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%   DATA:               Data structure; A data structure that contains at 
%                       least the fields specified by WHICHSIGFIELD,
%                       WHICHBAQFIELD. and WHICHFS. For example,
%                       465nm stream as data.sig, 405nm stream as data.baq,
%                       and sampling rate as data.fs.
% 
%   WHICHSIGFIELD:      String; The name of the field containing the
%                       signal stream. For example, 'sig'.
%
%   WHICHBAQFIELD:      String; The name of the field containing the
%                       background stream. For example, 'baq'.
%
%   WHICHFS:             String; The name of the field containing the 
%                       sampling rate of the streams. For example, 'fs'.
%
% OPTIONAL INPUTS:
%   BAQSCALINGTYPE:      String; Specifies the type of background scaling 
%                        to apply.
%                       'frequency': Scales the background to the signal 
%                           channel based on ratio of specified frequency  
%                           bands in the FFT of the channels.
%                       'sigmean': Scales the background to the signal 
%                           channel based on the ratio of the mean of the
%                           signal to the mean of the background.
%                       'OLS': Uses ordinary least-squares regression to
%                           generate scaled background.
%                       'detrendedOLS': Removes the linear trend from signal
%                           and background streams prior to using ordinary 
%                           least-squares regression to generate scaled
%                           background.
%                       'smoothedOLS': Applies lowess smoothing to the signal
%                           and background streams prior to using ordinary 
%                           least-squares regression to generate scaled
%                           background.
%                       'IRLS': Uses iteratively reweighted least squares
%                           regression to generate scaled background.
%                       Default: 'frequency'.
%
%   BAQSCALINGFREQMIN:  Numeric; Only used with 'frequency' baqscaling. 
%                       Frequency (Hz) minimum threshold for scaling the 
%                       background to signal channel. Frequencies above 
%                       this value will be included in the scaling factor 
%                       determination.
%                       Default: 10 Hz.
%
%   BAQSCALINGFREQMAX:  Numeric; Only used with 'frequency' baqscaling. 
%                       Frequency (Hz) maximum threshold for scaling the 
%                       background to signal channel. Frequencies below 
%                       this value will be included in the scaling factor 
%                       determination.
%                       Default: 100 Hz.
%
%   BAQSCALINGPERC:     Numeric; Only used with 'frequency' and 'sigmean' 
%                       baqscaling. Adjusts the background scaling factor 
%                       to be a percent of the derived scaling factor value. 
%                       Default: 1 (100%)
%
%   SUBTRACTIONOUTPUT:  String; Output type for the subtracted data.
%                       'dff':  delta F/F; ((sig/baq)./baq)
%                       'df':   delta F; (sig/baq).
%                       Default: 'dff'
%
%   ARTIFACTREMOVAL:    Numeric; Set to 1 to automatically detect and
%                       remove artifacts from data streams. Default: 0.
%
%   FILTERTYPE:         String; Specifies the type of filter to be applied
%                       to the subtracted data stream.
%                       'nofilter': No filter will be applied.
%                       'bandpass': A bandpass filter will be applied. 
%                       'highpass': Only the high pass filter will be 
%                                   applied.
%                       'lowpass':  Only the low pass filter will be 
%                                   applied.
%                       Default: 'bandpass'
%
%   PADDING:            Numeric; Defaults to 1, which applies padding. 
%                       Padding takes the specific percent of the start of
%                       the stream, flips it, and appends it to the data
%                       before filtering. Appended data is trimmed after 
%                       filtration. Set to 0 to turn off padding.
%                       Default: 1
%
%   PADDINGPERC:        Numeric; Percent of data length to use to determine the
%                       number of samples to be appended to the beginning
%                       and end of data in padding. Set to minimum 10%.
%                       Default: 0.1 (10%)
%   
%   FILTERORDER:        Numeric; The order to be used for the Butterworth
%                       filter. Default: 3.
%
%   HIGHPASSCUTOFF:     Numeric; The cutoff frequency (Hz) to be used for 
%                       the high pass Butterworth filter. Default: 2.2860.
%
%   LOWPASSCUTOFF:      Numeric; The cutoff frequency (Hz) to be used for 
%                       the low pass Butterworth filter. Default: 0.0051.
%
%   SUPRESSDISP:        Numeric; If set to anything other than 0, this will 
%                       suppress the command window displays. Default: 0
%
% OUTPUTS:
%       DATA:           Data structure; This is the original data structure 
%                       with added fields with the function inputs, scaled 
%                       background (baqscaled), subtracted signal (sigsub), 
%                       and subtracted and filtered signal (sigfilt).
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Settings
% Import optional inputs into a structure
    inputs = struct(...
        'whichsigfield',[],...
        'whichbaqfield',[],...
        'whichfs',[],...
        'baqscalingtype',[],...
        'baqscalingfreqmin',[],...
        'baqscalingfreqmax',[],...
        'baqscalingperc', [],...
        'subtractionoutput',[],...
        'artifactremoval',[],...
        'filtertype',[],...
        'padding', [],...
        'paddingperc', [],...
        'filterorder', [],...
        'highpasscutoff', [],...
        'lowpasscutoff',[],...
        'suppressdisp',[]);
    inputs = parseArgsLite(varargin,inputs);
    
    % Required inputs
    inputs.whichsigfield = whichsigfield;
    inputs.whichbaqfield = whichbaqfield;
    inputs.whichfs = whichfs;

    % Prepare defaults and check for optional inputs
    if isempty(inputs.baqscalingtype)
        baqscalingtype = 'frequency'; % Subtraction type defaults to frequency scaling
        inputs.baqscalingtype = baqscalingtype;
    else
        baqscalingtype = inputs.baqscalingtype;
    end
    if isempty(inputs.baqscalingfreqmin)
        baqscalingfreqmin = 10; % Background scaling frequency min defaults to 10 Hz
        inputs.baqscalingfreqmin = baqscalingfreqmin;
    else
        baqscalingfreqmin = inputs.baqscalingfreqmin;
    end
    if isempty(inputs.baqscalingfreqmax)
        baqscalingfreqmax = 100; % Background scaling frequency max defaults to 100 Hz
        inputs.baqscalingfreqmax = baqscalingfreqmax;
    else
        baqscalingfreqmax = inputs.baqscalingfreqmax;
    end
    if isempty(inputs.baqscalingperc)
        baqscalingperc = 1; % Background scaling percent defaults to 100%
        inputs.baqscalingperc = baqscalingperc;
    else
        baqscalingperc = inputs.baqscalingperc;
    end
    if isempty(inputs.subtractionoutput)
        subtractionoutput = 'dff'; % Subtraction output defaults to delta F/F
        inputs.subtractionoutput = subtractionoutput;
    else
        subtractionoutput = inputs.subtractionoutput;
    end
    if isempty(inputs.artifactremoval)
        artifactremoval = 0; % Subtraction output defaults to delta F/F
        inputs.artifactremoval = artifactremoval;
    else
        artifactremoval = inputs.artifactremoval;
    end
    if isempty(inputs.filtertype)
        filtertype = 'bandpass'; % Filter type defaults to bandpass filter
        inputs.filtertype = filtertype;
    else
        filtertype = inputs.filtertype;
    end
    if isempty(inputs.padding)
        padding = 1; % Pre-filter padding defaults to apply padding (1)
        inputs.padding = padding;
    else
        padding = inputs.padding;
    end
    if isempty(inputs.paddingperc)
        paddingperc = 0.1; % Pre-filter padding length defaults to 10% of stream length
        inputs.paddingperc = paddingperc;
    else
        paddingperc = inputs.paddingperc;
    end
    if isempty(inputs.filterorder)
        filterorder = 3; % Filter order defaults to 3rd order
        inputs.filterorder = filterorder;
    else
        filterorder = inputs.filterorder;
    end
    if isempty(inputs.highpasscutoff)
        highpasscutoff = .0051; % High pass frequency cuttoff defaults to 0.0051 Hz
        inputs.highpasscutoff = highpasscutoff;
    else
        highpasscutoff = inputs.highpasscutoff;
    end
    if isempty(inputs.lowpasscutoff)
        lowpasscutoff = 2.2860; % Low pass frequency cuttoff defaults to 2.2860 Hz
        inputs.lowpasscutoff = lowpasscutoff;
    else
        lowpasscutoff = inputs.lowpasscutoff;
    end
    if isempty(inputs.suppressdisp)
        suppressdisp = 0; % Display text defaults to showing
        inputs.suppressdisp = 0;
    else
        suppressdisp = inputs.suppressdisp;
    end

%% Display subtraction settings
    if suppressdisp == 0
        if strcmp(baqscalingtype,'frequency')==true % Display for frequency background scaling - default
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with frequency domain scaling (threshold at ', num2str(baqscalingfreqmin), 'hz). Subtracted signal will be output as ', subtractionoutput, ' to data.sigsub.'));
        elseif strcmp(baqscalingtype,'sigmean')==true % Display for signal mean background scaling
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with signal mean scaling. Subtracted signal will be output as ', subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(baqscalingtype,'OLS')==true % Display for ordinary least squares regression background scaling
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with ordinary least-squares regression (OLS). Subtracted signal will be output as ', subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(baqscalingtype,'detrendedOLS')==true % Display for linear detrending and ordinary least squares regression background scaling
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with detrending and ordinary least-squares regression (OLS). Subtracted signal will be output as ', subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(baqscalingtype,'smoothedOLS')==true % Display for time domain subtraction
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with lowess smoothing and ordinary least-squares regression (OLS). Subtracted signal will be output as ', subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        elseif strcmp(baqscalingtype,'IRLS')==true % Display for time domain subtraction
            disp(append('SUBTRACTFPDATA: ','Subtracting ',whichbaqfield, ' from ',whichsigfield, ' with iteratively reweighted least squares regression (IRLS). Subtracted signal will be output as ', subtractionoutput, ' to data.sigsub.'));
            disp(append('WARNING: ','Deviation from PASTa protocol default in baqscalingtype.'));
        else
            disp(append('ERROR: Baq scaling type issue - baqscalingtype set to ', baqscalingtype, '. Function inputs limited to frequency, sigmean, OLD, detrendedOLS, smoothedOLS, or IRLS.'));
        end

        if strcmp(filtertype,'nofilter')==true % Display for filter settings
            disp(append('WARNING: filtertype manually set to nofilter - NO FILTER APPLIED'));
        else
            disp(append('Filter type set to ', filtertype, '. Subtracted and filtered signal will be output as data.sigfilt.'));
            if padding == 1 % Padding application display
                disp('   NOTE: Padding applied prior to filtering and removed from final output.')
            end
        end
        if artifactremoval == 1 % Display if artifact removal is turned on
            disp(append('Artifact removal set to 1 - artifacts will be automatically detected and replaced with NaNs in subtracted and filtered outputs.'));
        end

        disp('INPUTS:') % Display all input values
        disp(inputs)
    end

%% Subtract and Filter Data
    for eachfile = 1:length(data)
        fs = data(eachfile).(whichfs);

        % Prepare filters
        if strcmp(filtertype,'bandpass')==true
            highpasscutoffval = round(highpasscutoff/floor(fs/2),6);
            lowpasscutoffval = round(lowpasscutoff/floor(fs/2),6);
            [a,b] = butter(filterorder,highpasscutoffval,'high');
            [c,d] = butter(filterorder,lowpasscutoffval,'low');
        elseif strcmp(filtertype, 'highpass')==true
            highpasscutoffval = round(highpasscutoff/floor(fs/2),6);
            [a,b] = butter(filterorder,highpasscutoffval,'high');
        elseif strcmp(filtertype, 'lowpass')==true
            lowpasscutoffval = round(lowpasscutoff/floor(fs/2),6);
            [c,d] = butter(filterorder,lowpasscutoffval,'low');
        end

        data(eachfile).inputs_subtractFPdata = inputs; % Add inputs to data frame
        if suppressdisp == 0
            fprintf('Subtracting file number: %.f \n',eachfile) % Display which file is being subtracted
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
            if strcmp('frequency',baqscalingtype)==true % Frequency domain scaling - scales to power
                baq_centered = baq - mean(data(eachfile).(whichbaqfield)); % Center background at 0
                sig_centered = sig - mean(data(eachfile).(whichsigfield)); % Center signal at 0

                [sigFFT, sigF] = preparestreamFFT(sig_centered,fs); % Prep FFT
                [baqFFT, baqF] = preparestreamFFT(baq_centered,fs); % Prep FFT

                sigFidxs = sigF > baqscalingfreqmin & sigF < baqscalingfreqmax; % Find indices of frequencies above the set min threshold and below the set max threshold
                baqFidxs = baqF > baqscalingfreqmin & baqF < baqscalingfreqmax; % Find indices of frequencies above the set min threshold and below the set max threshold

                sigPower = sigFFT(sigFidxs).^2; % Compute average power in the selected band
                baqPower = baqFFT(baqFidxs).^2; % Compute average power in the selected band

                powerRatio = mean(sigPower) / mean(baqPower); % Find the ratio of mean powers between signal and background (same as RMS)
                baqscalingfactor = sqrt(powerRatio) * baqscalingperc;  % or remove baqscalingperc if not needed

                baqscaled = (baq_centered*baqscalingfactor) + mean(sig); % Adjust back to same units as raw signal
                data(eachfile).baqscalingfactor = baqscalingfactor; % Add constant baqscalingfactor to data

                if baqscalingfactor > 3 % Display a warning if the baqscaling factor is high
                    disp(append('   WARNING: Possible over scaling. baqscalingfactor = ',num2str(baqscalingfactor)))
                end
            elseif strcmp('frequency_amplitude',baqscalingtype)==true % Frequency domain scaling - scales to amplitude
                baq_centered = baq - mean(data(eachfile).(whichbaqfield)); % Center background at 0
                sig_centered = sig - mean(data(eachfile).(whichsigfield)); % Center signal at 0

                [sigFFT, sigF] = preparestreamFFT(sig_centered,fs); % Prep FFT
                [baqFFT, baqF] = preparestreamFFT(baq_centered,fs); % Prep FFT

                sigFidxs = sigF > baqscalingfreqmin; % Find indices of frequencies above the set threshold
                baqFidxs = baqF > baqscalingfreqmin; % Find indices of frequencies above the set threshold

                baqscalingfactor = (mean(sigFFT(sigFidxs))/mean(baqFFT(baqFidxs)))*baqscalingperc; % Find power ratio of signal to background for frequencies above the set threshold
                baqscaled = (baq_centered*baqscalingfactor) + mean(sig); % Adjust back to same units as raw signal

                data(eachfile).baqscalingfactor = baqscalingfactor; % Add constant baqscalingfactor to data
            elseif strcmp('sigmean',baqscalingtype)==true % Signal mean time domain scaling
                baqscalingfactor = (mean(sig)/(mean(baq)))*baqscalingperc; % Find the scaling factor for background to signal based on the ratio of the means
                baqscaled = baq*baqscalingfactor;
                data(eachfile).baqscalingfactor = baqscalingfactor; % Add constant baqscalingfactor to data
            elseif strcmp('OLS',baqscalingtype)==true % OLS regression scaling; Per GuPPY, see Sherathiya et al 2021, https://www.nature.com/articles/s41598-021-03626-9
                bls=polyfit(baq,sig,1);
                baqscaled=bls(1).*baq+bls(2);      
            elseif strcmp('detrendedOLS',baqscalingtype)==true % Linear detrend and OLS regression scaling
                sigmean = mean(sig);
                baq = detrend(baq); % Detrend background
                sig = detrend(sig)+sigmean; % Detrend signal
                bls = polyfit(baq,sig,1); % Fit background to signal with least-squares linear regression
                baqscaled = (bls(1).*baq)+ bls(2); % Fit to signal
                data(eachfile).baqdetrend = baq; % Add detrended signal to data
                data(eachfile).sigdetrend = sig; % Add detrended background to data
            elseif strcmp('smoothedOLS',baqscalingtype)==true % Lowess smoothing and OLS regression scaling; Per pMAT, see Bruno et al 2021, https://www.sciencedirect.com/science/article/abs/pii/S0091305720307413?via=ihub#f0020
                baq=smooth(baq,0.002,'lowess')'; 
                sig=smooth(sig,0.002,'lowess')';
                bls=polyfit(baq(1:end),sig(1:end),1);
                baqscaled=bls(1).*baq+bls(2);
                data(eachfile).baqsmoothed = baq; % Add detrended signal to data
                data(eachfile).sigsmoothed = sig; % Add detrended background to data
            elseif strcmp('IRLS',baqscalingtype)==true % IRLS regression scaling; See Keevers et al 2024, https://www.researchsquare.com/article/rs-3549461/v2
                IRLS_coeffs = reshape(flipud(robustfit(baq, sig, 'bisquare', 1.4, 'on')), [1, 2]);
                baqscaled = polyval(IRLS_coeffs,baq);
            else
                disp(append('ERROR: Baq scaling type issue - baqscalingtype set to ', baqscalingtype, '. Function inputs limited to frequency, sigmean, OLS, detrendedOLS, smoothedOLS, or IRLS.'));
            end
            data(eachfile).(append(whichbaqfield,'scaled')) =  baqscaled; % Add scaled background to data frame

        %% Subtract data
            if strcmp(subtractionoutput,'dff')==true % Output delta F/F
                data(eachfile).sigsub = (sig-baqscaled)./baqscaled*100;
            elseif strcmp(subtractionputput,'df')==true % Output delta F
                data(eachfile).sigsub = (sig-baqscaled);
            else
                disp(append('ERROR: Subtraction Output type issue - subtractionoutput set to ', subtractionoutput, '. Function inputs limited to df or dff.'));
            end
           
        %% OPTIONAL: Remove Artifacts
            if artifactremoval == 1
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
            if strcmp(filtertype, 'nofilter')==true % Skip filtering
            else
                if padding==1 % Filter with padding
                    nsamplesedge = floor(length(data(eachfile).sigsub)*paddingperc); % Padding: determine number of samples to append to start and end
                    firstsamples = [flip(data(eachfile).sigsub(1:nsamplesedge),2)]; % Padding: Extract data to append to beginning of signal
                    lastsamples = [flip(data(eachfile).sigsub(length(data(eachfile).sigsub)-(nsamplesedge-1):length(data(eachfile).sigsub)),2)]; % Padding: Extract data to append to end of signal
                    sigsubdata = [firstsamples data(eachfile).sigsub lastsamples]; % Padding: Create padded data
                    
                    % Apply filter
                    if strcmp(filtertype,'bandpass') 
                        filtdata = filtfilt(a,b,sigsubdata); % Apply the highpass
                        filtdata = filtfilt(c,d,filtdata); % Apply the lowpass
                    elseif strcmp(filtertype,'highpass')
                        filtdata = filtfilt(a,b,sigsubdata); % Apply the highpass
                    elseif strcmp(filtertype,'lowpass')
                        filtdata = filtfilt(c,d,sigsubdata); % Apply the highpass
                    else
                        disp('WARNING: Filter type undefined - no filter applied.')
                    end
                    data(eachfile).sigfilt = filtdata(nsamplesedge+1:length(filtdata)-nsamplesedge); % Padding: Remove padding from beginning and end
                else % Filter without padding
                    if strcmp(filtertype,'bandpass')
                        filtdata = filtfilt(a,b,data(eachfile).sigsub); % Apply the highpass
                        filtdata = filtfilt(c,d,filtdata); % Apply the lowpass
                    elseif strcmp(filtertype,'highpass')
                        filtdata = filtfilt(a,b,data(eachfile).sigsub); % Apply the highpass
                    elseif strcmp(filtertype,'lowpass')
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
            if artifactremoval == 1
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