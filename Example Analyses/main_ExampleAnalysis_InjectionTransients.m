%% EXAMPLE ANALYSIS
% PROJECT SUMMARY: This is an example analysis to determine changes in GCaMP6f transients in VTA dopamine neurons after within session saline 
% or morphine injection. Each subject has two recording sessions: a saline control session and a morphine administration session. Each session
% recording consists of a 15 minute pre-injection baseline, injection, and a 60 minute post injection period.


%%%%%%%%%%%%%%%%%%%%%%%%%%%% TABLE OF CONTENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE OF CONTENTS
% #1  Basic Protocol 2.3 - Prepare the MATLAB environment
% #2  Basic Protocol 2.4 - Create experimentkey
% #3  Basic Protocol 2.5 - Extract data
% #4  Basic Protocol 2.6 - Load fiber photometry data into structure
% #5  Basic Protocol 2.7 - Prepare session start and end indices
% #6  Basic Protocol 2.8 - Crop fiber photometry data streams
% #7  Basic Protocol 3.1 - Subtract and filter fiber photometry data
% #8  Basic Protocol 3.2 - Plot stream traces
% #9  Basic Protocol 3.4 - Plot stream frequency spectra
% #10 Basic Protocol 3.6 - Normalize subtracted and filtered data stream
% #11 Basic Protocol 3.7 - Plot normalized streams
% #12 Basic Protocol 4.1 - Determine transient detection thresholds
% #13 Basic Protocol 4.2 - Find and quantify transient events
% #14 Basic Protocol 4.3 - Bin transient events
% #15 Basic Protocol 4.4 - Export transient events
% #16 Basic Protocol 4.5 - Plot whole session transients
% #17 Basic Protocol 4.6 - Plot binned transients
% #18 Basic Protocol 4.7 - Plot whole session transient traces
% #19 Basic Protocol 4.8 - Plot binned transient traces
% #20 Basic Protocol 4.9 - Summarize transients by session
% #21 Basic Protocol 4.10 - Summarize transients by bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% #1  Basic Protocol 2.3 - Prepare the MATLAB environment
% Set up user path inputs
rootdirectory = getenv("USERPROFILE") + "\";
rootdirectory_files = 'C:\Users\rmdon\Box\PASTaExampleFiles\Example Analyses\'; % Path for analysis files - this is where the data and keys are saved

analysisfolder = 'C:\Users\rmdon\Box\PASTaExampleFiles\Example Analyses\Injection Transients\Analysis\'; % Folder to output analysis csv files to
figurefolder = 'C:\Users\rmdon\Box\PASTaExampleFiles\Example Analyses\Injection Transients\Figures\'; % Folder to output figures to

% Add data folders to MATLAB path
addpath(genpath(rootdirectory_files));

% Add PASTa repository folders to MATLAB path - only needed if PASTa is not installed as a MATLAB Toolbox 
% through the MATLAB Toolbox Add-On Explorer or the .mbtlx file in the PASTa repository.
addpath(genpath('C:\Users\rmdon\OneDrive\Desktop\GitHub_Repositories\PASTa\')); % Set this to the path to the folder location of the local copy of PASTa functions

% Load in experiment key names - Subject Key and File Key
subjectkeyname = 'SubjectKey_ExampleAnalysis_MorphineTransients.csv'; % Name of csv file containing subject information; set to '' if not using a Subject Key
filekeyname = 'FileKey_ExampleAnalysis_MorphineTransients.csv'; % Name of csv file containing session information and paths

%% #2  Basic Protocol 2.4 - Create experimentkey
% Load subject key and file key into a data structure and append rootdirectory to RawFolderPath and ExtractedFolderPaths
[experimentkey] = createExperimentKey(rootdirectory_files, subjectkeyname, filekeyname);

%% #3  Basic Protocol 2.5 - Extract data
% Extract raw data from blocks using the function 'extractTDTdata' to extract raw data blocks.
% Set up inputs
sigstreamnames = {'x65A', '465A'}; % All names of signal streams across files
baqstreamnames = {'x05A', '405A'}; % All names of background streams across files
rawfolderpaths = string({experimentkey.RawFolderPath})'; % Create string array of raw folder paths
extractedfolderpaths = string({experimentkey.ExtractedFolderPath})'; % Create string array of extracted folder paths

extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames); % extract data

%% #4  Basic Protocol 2.6 - Load fiber photometry data into structure
% Load previously extracted data blocks and tie to experiment key. 
% Each block is loaded as a row in the data structure.
[rawdata] = loadKeydata(experimentkey); % Load data based on the experiment key into the structure 'rawdata'
 
%% #5  Basic Protocol 2.7 - Prepare session start and end indices
% Prepare session start and end indices
preinjectionlength = 15; % Minutes pre injection
postinjectionlength = 60; % Minutes post injection
for eachfile = 1:length(rawdata)
    rawdata(eachfile).sessionstart = rawdata(eachfile).injt(1) - (floor(preinjectionlength*60*rawdata(eachfile).fs)); % Find and save start index to rawdata structure field
    sessionend = rawdata(eachfile).injt(2) + (floor(postinjectionlength*60*rawdata(eachfile).fs)); % Find end index and check that it's within the length of the full session
    if length(rawdata(eachfile).sig) >= sessionend % Save end index to rawdata structure field
        rawdata(eachfile).sessionend = sessionend;
    else  % If the found end index is greater than the length of the signal, use the length of the signal
        disp(append('WARNING: Session short for file ', num2str(eachfile))) 
        rawdata(eachfile).sessionend = length(rawdata(eachfile).sig);
    end
end

%% #6  Basic Protocol 2.8 - Crop fiber photometry data streams
% Crop data: remove pre and post experimental session samples.
% NOTE: Cropping is the only part of the pipeline that will alter the loaded data fields. 
% To ensure you only complete this step once per analysis, it is reccomended to input structure 
% 'rawdata' to the function and output a new structure 'data'.
cropstartfieldname = 'sessionstart'; % name of field with session start index
cropendfieldname = 'sessionend'; % name of field with session end index
streamfieldnames = {'sig', 'baq'}; % which streams to crop
epocsfieldnames = {'injt','sess'}; % which epocs to adjust to maintain relative position - OPTIONAL INPUT

[data] = cropFPdata(rawdata,cropstartfieldname,cropendfieldname, streamfieldnames,'epocsfieldnames', epocsfieldnames); % Output cropped data into new structure called data

%% #7  Basic Protocol 3.1 - Subtract and filter fiber photometry data
% Subtract and filter data with default settings
sigfieldname = 'sig';
baqfieldname = 'baq';
fsfieldname = 'fs';

[data] = subtractFPdata(data,sigfieldname,baqfieldname,fsfieldname); % adds sigsub (subtracted stream) and sigfilt (subtracted and filtered stream) to data frame

%% #8  Basic Protocol 3.2 - Plot stream traces
% Plot whole session streams for each file. 
% Use plotTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
% Manually save plots to allow for customization (addition of injection start/stop lines)
outputfiletype = 'png';

for eachfile = 1:length(data)
    fileindex = eachfile;
    maintitle = append(num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'SessionTraces_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);

    alltraces = plotTraces(data,fileindex,maintitle); % Save plot into object for customization
    for eachtile = 1:5 % Add injection start/stop lines to each stream tile
        nexttile(eachtile)
        xline(data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',8)
        xline(data(eachfile).injt(2),'--','Color','#C40300','FontSize',8)
    end    

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]); % Manually save the figure
    exportgraphics(gcf,append(plotfilepath,'.',outputfiletype),'Resolution',300)
end

%% #9  Basic Protocol 3.4 - Plot stream frequency spectra
% Plot whole session FFT power plots for each file
% Use plotFFTpower to plot all frequency magnitude plots - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
fsfieldname = 'fs'; % Prepare field names for function inputs 

for eachfile = 1:length(data) % Plot each file
    fileindex = eachfile;
    maintitle = append(num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'SessionFFTpower_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);

    plotFFTpower(data,fileindex,maintitle,fsfieldname,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath); % Autmatically save plots to specified file path
end

% Plot whole session FFT magnitude plots for each file
% Use plotFFTmag to plot all frequency magnitude plots - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
for eachfile = 1:length(data) % Plot each file
    fileindex = eachfile;
    maintitle = append(num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'SessionFFTmag_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);

    plotFFTmag(data,fileindex,maintitle,fsfieldname,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
end

%% #10 Basic Protocol 3.6 - Normalize subtracted and filtered data stream
% Option 1: Normalize to whole session mean and standard deviation
streamfieldname = 'sigfilt'; % Prepare field names for function inputs 
[data] = normSession(data,streamfieldname); % Outputs Z scored stream based on whole session mean and SD

% Option 2: Normalize to session baseline mean and standard deviation
for eachfile = 1:length(data) % Prepare indexes for baseline period start and end
    data(eachfile).BLstart = 1;
    data(eachfile).BLend =  data(eachfile).injt(1);
end

BLstartfieldname = 'BLstart'; % Prepare field names for function inputs 
BLendfieldname = 'BLend';
[data] = normBaseline(data,streamfieldname,BLstartfieldname,BLendfieldname);  % Outputs Z scored stream based on session baseline (pre-injection) mean and SD

%% EXPERIMENT SPECIFIC: Remove injection time window from normalized signal
% This loop removes samples between the start and end of the injection. For
% cropped data, injt(1) will be used as the injection time point.
normstreams = {'sigfilt', 'sigfiltz_normsession','sigfiltz_normbaseline'}; % List of all streams to remove the injection period from
for eachfile = 1:length(data)
    for eachstream = 1:length(normstreams)
        currstream = char(normstreams(eachstream)); % Pull stream name from cell array into a character object

        allindices = (1:length(data(eachfile).(currstream))); % Temporary list of all indices in the stream
        includeindices = (allindices < data(eachfile).injt(1) | allindices > data(eachfile).injt(2)); % Find indices before injection start and after injection end
        data(eachfile).(append(currstream,'_injcropped')) = data(eachfile).(currstream)(includeindices); % Add new fields ending in '_trimmed' to data structure
    end 
end

%% #11 Basic Protocol 3.7 - Plot normalized streams
% Use plotNormTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
streams = {'sigfiltz_normsession', 'sigfiltz_normbaseline'};
streamtitles = {'Normalized to Whole Session', 'Normalized to Baseline'};
outputfiletype = 'png';

for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'SessionNormTraces_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);
    normtraces = plotNormTraces(data,eachfile,streams,'fs',maintitle,streamtitles);


    for eachtile = 1:2
        nexttile(eachtile)
        xline(data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',8)
    end    

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]);
    exportgraphics(gcf,append(plotfilepath,'.',outputfiletype),'Resolution',300)
end


%% #12 Basic Protocol 4.1 - Determine transient detection thresholds
% Prepare thresholds for both normalized and df/f streams
for eachfile = 1:length(data)
    % Add thresholds
    data(eachfile).SDthreshold = 2.6;
    data(eachfile).SDthresholdsigfilt = 2.6*std(data(eachfile).sigfilt_injcropped(1:data(eachfile).injt(1)));
    
    % Add mean and SD of subtracted and filtered signal to data structure in case of alternative threshold setting needs
    data(eachfile).mean_sigfilt = mean(data(eachfile).sigfilt_injcropped(1:data(eachfile).injt(1)));
    data(eachfile).SD_sigfilt = std(data(eachfile).sigfilt_injcropped(1:data(eachfile).injt(1))); 
end

%% #13 Basic Protocol 4.2 - Find and quantify transient events
% Create list of variables to add to the new data structure with the output transient events
addvariablesfieldnames = [fieldnames(experimentkey); {'params'}]; % This makes a list of all fieldnames in the experimentkey and adds the 'params' field

% Find transients based on pre-peak baseline window mean - reccomended as the first pass choice for transient analysis
[transientdata_sigfiltznormBL] = findTransients(data,addvariablesfieldnames,'sigfiltz_normbaseline_injcropped','SDthreshold','fs');
[transientdata_sigfilt] = findTransients(data,addvariablesfieldnames,'sigfilt_injcropped','SDthresholdsigfilt','fs');


%% #14 Basic Protocol 4.3 - Bin transient events
% Bin transients into time bins - default 5 mins
[transientdata_sigfiltznormBL] = binTransients(transientdata_sigfiltznormBL);
[transientdata_sigfilt] = binTransients(transientdata_sigfilt);

% Bin transients with 3 minute bins
[transientdata_sigfiltznormBL] = binTransients(transientdata_sigfiltznormBL,'binlengthmins',3);
[transientdata_sigfilt] = binTransients(transientdata_sigfilt,'binlengthmins',3);

% OPTIONAL: Create custom bins
% If needed, transients can be binned into varying bins based on trial start and end indexes, or other relevant time points. 
% Start and end indexes for each desired bin must be prepared and then passed to the function binTransients as optional inputs. 
% This is an example of custom bin start and end index creation based on time points - in this case, just pre and post injection.
for eachfile = 1:length(data)
    % Prepare pre-injection start and end indexes
    transientdata_sigfiltznormBL(eachfile).binstart(1) = 1;
    transientdata_sigfiltznormBL(eachfile).binend(1) = data(eachfile).injt(1);

    transientdata_sigfilt(eachfile).binstart(1) = 1;
    transientdata_sigfilt(eachfile).binend(1) = data(eachfile).injt(1);

    % Prepare post-injection start and end indexes
    transientdata_sigfiltznormBL(eachfile).binstart(2) = transientdata_sigfiltznormBL(eachfile).binend(1)+1;
    transientdata_sigfiltznormBL(eachfile).binend(2) = length(data(eachfile).sigfiltz_normbaseline_injcropped);

    transientdata_sigfilt(eachfile).binstart(2) = transientdata_sigfilt(eachfile).binend(1)+1;
    transientdata_sigfilt(eachfile).binend(2) = length(data(eachfile).sigfilt_injcropped);
end

% Prep field names in variables
binstartfieldname = 'binstart';
binendfieldname = 'binend';

% Bin transients with custom bin size
[transientdata_sigfiltznormBL] = binTransients(transientdata_sigfiltznormBL,'manuallydefinebins',1,'binstartfieldname',binstartfieldname,'binendfieldname',binendfieldname);
[transientdata_sigfilt] = binTransients(transientdata_sigfilt,'manuallydefinebins',1,'binstartfieldname',binstartfieldname,'binendfieldname',binendfieldname);

%% #15 Basic Protocol 4.4 - Export transient events
% Saves all individual transient events to one table and exports the table to a csv file
addvariables = {'SubjectID','TreatNum','InjType','Weight','FiberPlacement','Virus'};
alltransients_Z = exportTransients(transientdata_sigfiltznormBL,'transientquantification',analysisfolder,addvariables,'exportfilename','transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold.csv');


%% #16 Basic Protocol 4.5 - Plot whole session transients
% Use plotTransients to plot the whole session trace with detected transients for each file.
for eachfile = 1:length(data)
    fileindex = eachfile;
    maintitle = append('Subject ',num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'SessionTransients_blmin_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);

    sessiontransients = plotTransients(data,fileindex,'sigfiltz_normsession_injcropped','fs',transientdata_sigfiltznormBL,maintitle,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
end

%% #17 Basic Protocol 4.6 - Plot binned transients
% Use plotTransientBins to plot session bin traces with detected transients for each file.
for eachfile = 1:length(data)
    fileindex = eachfile;
    maintitle = append('Subject ',num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'SessionBins_blmin_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);

    allbins = plotTransientBins(data,fileindex,'sigfiltz_normsession_injcropped',transientdata_sigfiltznormBL,'Bin_5min',maintitle,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
end

%% #18 Basic Protocol 4.7 - Plot whole session transient traces
% Use plotTransientTraces to plot overlaid traces of all detected transients in the session for each file.
for eachfile = 1:length(data)
    fileindex = eachfile;
    maintitle = append('Subject ',num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'TransientTraces_blmin_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);
    alltransienttraces = plotTransientTraces(transientdata_sigfiltznormBL,fileindex,maintitle,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
end


%% #19 Basic Protocol 4.8 - Plot binned transient traces
% Use plotTransientTraces to plot overlaid traces of detected transients by bin for each file.
for eachfile = 1:length(data)
    fileindex = eachfile;
    maintitle = append('Subject ',num2str(data(eachfile).SubjectID),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    plotfilepath = append(figurefolder,'TransientTracesby5minBin_blmin_',num2str(data(eachfile).SubjectID),'_',data(eachfile).InjType);

    alltransienttracebins = plotTransientTraceBins(transientdata_sigfiltznormBL,fileindex,'Bin_5min',maintitle,'saveoutput',1,'outputfiletype','png','plotfilepath',plotfilepath);
end


%% #20 Basic Protocol 4.9 - Summarize transients by session
% summarizeTransients finds the whole session frequency (in total peaks, peaks per minute, and peak per second/hz), mean values for quantification
% variables (amp, rise, fall, width, AUC, etc), and total number of compound events.
[transientdata_sigfiltznormBL] = summarizeTransients(transientdata_sigfiltznormBL);

%summarizeBinTransients finds the per bin values within each session for frequency (in total peaks, peaks per minute, and peak per second/hz), 
% mean values for quantification variables (amp, rise, fall, width, AUC, etc), and total number of compound events.
[transientdata_sigfiltznormBL] = summarizeBinTransients(transientdata_sigfiltznormBL,'Bin_5min');


%% #21 Basic Protocol 4.10 - Summarize transients by bin
% summarizeBinTransients finds the per bin values within each session for frequency (in total peaks, peaks per minute, and peak per second/hz), 
% mean values for quantification variables (amp, rise, fall, width, AUC, etc), and total number of compound events.
[transientdata_sigfiltznormBL] = summarizeBinTransients(transientdata_sigfiltznormBL,'Bin_5min');


% Export transients with added fields for subject and treatment using the
% EXPORTSESSIONTRANSIENTS function. See Basic Protocol 4.4.
addvariables = {'SubjectID','TreatNum','InjType','Weight','FiberPlacement','Virus'};
alltransients_Z_session = exportTransients(transientdata_sigfiltznormBL,'transientsummary_session',analysisfolder,addvariables,'exportfilename','transientsummary_session_sigfiltz_normbaseline_injcropped_SDthreshold.csv');
alltransients_Z_bin = exportTransients(transientdata_sigfiltznormBL,'transientsummary_Bin_5min',analysisfolder,addvariables,'exportfilename','transientsummary_Bin5min_sigfiltz_normbaseline_injcropped_SDthreshold.csv');
