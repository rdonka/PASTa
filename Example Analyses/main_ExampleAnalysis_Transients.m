%% EXAMPLE ANALYSIS
% PROJECT SUMMARY: This is an example analysis to determine changes in GCaMP6f transients in VTA dopamine neurons after within session saline 
% or morphine injection. Each subject has two recording sessions: a saline control session and a morphine administration session. Each session
% recording consists of a 15 minute pre-injection baseline, injection, and a 60 minute post injection period.

%% Set up paths and analysis keys
% Set up user path inputs
computeruserpath =  'C:\Users\rmdon\'; % Computer user unique portion of file path
analysisfolder = 'Box\PASTaExampleFiles\Example Analyses\Injection Transients\Analysis\'; % Folder to output analysis csv files to
figurefolder = 'Box\PASTaExampleFiles\Example Analyses\Injection Transients\Figures\'; % Folder to output figures to

% Create full paths with computeruserpath appended
analysispath = append(computeruserpath,analysisfolder); 
figurepath = append(computeruserpath,figurefolder);

% Add GitHub Repositories and data folders to MATLAB path
addpath(genpath(append(computeruserpath,'Onedrive\Desktop\GitHub_Repositories\PASTa\'))); % Path for GitHub repository
addpath(genpath(append(computeruserpath,'Box\PASTaExampleFiles\'))); % Path for analysis files - this is where the keys are saved\\

% Load in experiment key names - Subject Key and File Key
subjectkeyname = 'SubjectKey_ExampleAnalysis_MorphineTransients.csv'; % Name of csv file containing subject information; set to '' if not using a Subject Key
filekeyname = 'FileKey_ExampleAnalysis_MorphineTransients.csv'; % Name of csv file containing session information and paths

%% Load keys
% Load subject key and file key into a data structure and append computeruserpath to RawFolderPath and ExtractedFolderPaths
[experimentkey] = loadKeys(computeruserpath, subjectkeyname, filekeyname);

%% Extract data
% Extract raw data from blocks using the function 'extractTDTdata' to extract raw data blocks.
% Set up inputs
sigstreamnames = {'x65A', '465A'}; % All names of signal streams across files
baqstreamnames = {'x05A', '405A'}; % All names of background streams across files
rawfolderpaths = string({experimentkey.RawFolderPath})'; % Create string array of raw folder paths
extractedfolderpaths = string({experimentkey.ExtractedFolderPath})'; % Create string array of extracted folder paths

extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames); % extract data

%% Load data
% Load previously extracted data blocks and tie to experiment key. 
% Each block is loaded as a row in the data structure.
[rawdata] = loadKeydata(experimentkey); % Load data based on the experiment key into the structure 'rawdata'
 
%% Crop data
% Prepare session start and end indices
preinjectionlength = 15; % Minutes pre injection
postinjectionlength = 60; % Minutes post injection
for eachfile = 1:length(rawdata)
    sessionstart = rawdata(eachfile).injt(1) - (floor(preinjectionlength*60*rawdata(eachfile).fs)); % Find start index
    sessionend = rawdata(eachfile).injt(2) + (floor(postinjectionlength*60*rawdata(eachfile).fs)); % Find end index

    rawdata(eachfile).sessionstart = sessionstart; % Save start index to rawdata structure field

    if length(rawdata(eachfile).sig) >= sessionend % Save end index to rawdata structure field
        rawdata(eachfile).sessionend = sessionend;
    else  % If the found end index is greater than the length of the signal, use the length of the signal
        disp(append('WARNING: Session short for file ', num2str(eachfile))) 
        rawdata(eachfile).sessionend = length(rawdata(eachfile).sig);
    end
end

% Crop data: remove pre and post experimental session samples.
% NOTE: Cropping is the only part of the pipeline that will alter the loaded data fields. 
% To ensure you only complete this step once per analysis, it is reccomended to input structure 
% 'rawdata' to the function and output a new structure 'data'.
cropstart = 'sessionstart'; % name of field with session start index
cropend = 'sessionend'; % name of field with session end index
whichstreams = {'sig', 'baq'}; % which streams to crop
whichepocs = {'injt','sess'}; % which epocs to adjust to maintain relative position - OPTIONAL INPUT

[data] = cropFPdata(rawdata,cropstart,cropend, whichstreams,'whichepocs', whichepocs); % Output cropped data into new structure called data

%% Process data
% Subtract and filter data with default settings
sigfield = 'sig';
baqfield = 'baq';
fsfield = 'fs';

[data] = subtractFPdata(data,sigfield,baqfield,fsfield); % adds sigsub and sigfilt to data frame

%% Plot whole session streams for each file
% Use plotTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    alltraces = plotTraces(data,eachfile,maintitle);
    for eachtile = 1:5
        nexttile(eachtile)
        xline(data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',8)
        xline(data(eachfile).injt(2),'--','Color','#C40300','FontSize',8)
    end    

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]);
    plotfilepath = append(figurepath,'SessionTraces_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Normalize data
% To normalize to session mean:
[data] = normSession(data,'sigfilt'); % Outputs whole session z score

% To normalize to a session baseline:
for eachfile = 1:length(data) % prepare indexes for baseline period start and end
    data(eachfile).BLstart = 1;
    data(eachfile).BLend =  data(eachfile).injt(1);
end
[data] = normBaseline(data,'sigfilt','BLstart','BLend');

%% Remove injection time window from normalized signal
% This loop removes samples between the start and end of the injection. For
% cropped data, injt(1) will be used as the injection time point.
normstreams = {'sigfilt', 'sigfiltz_normsession','sigfiltz_normbaseline'};

for eachfile = 1:length(data)
    for eachstream = 1:length(normstreams)
        currstream = char(normstreams(eachstream)); % Pull stream name from cell array into a character object

        allindices = (1:length(data(eachfile).(currstream))); % Temporary list of all indices in the stream
        includeindices = (allindices < data(eachfile).injt(1) | allindices > data(eachfile).injt(2)); % Find indices before injection start and after injection end
        data(eachfile).(append(currstream,'_injcropped')) = data(eachfile).(currstream)(includeindices); % Add new fields ending in '_trimmed' to data structure
    end 
end

%% Find session transients
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients based on pre-peak baseline window minimum - reccomended as the first pass choice for transient analysis
[data] = findSessionTransients(data,'blmin','sigfiltz_normbaseline_injcropped','threshold3SD','fs');

% Find session transients based on pre-peak baseline window mean
[data] = findSessionTransients(data,'blmean','sigfiltz_normbaseline_injcropped','threshold3SD','fs');

% Find session transients based on pre-peak local minimum (last minumum before the peak in the baseline window)
[data] = findSessionTransients(data,'localmin','sigfiltz_normbaseline_injcropped','threshold3SD','fs');

% Bin session transients
[data] = binSessionTransients(data,'sigfiltz_normbaseline_injcropped','fs','sessiontransients_blmin_threshold3SD');
[data] = binSessionTransients(data,'sigfiltz_normbaseline_injcropped','fs','sessiontransients_blmean_threshold3SD');
[data] = binSessionTransients(data,'sigfiltz_normbaseline_injcropped','fs','sessiontransients_localmin_threshold3SD');


% Export transients with added fields for subject and treatment using the EXPORTSESSIONTRANSIENTS function
addvariables = {'Subject','TreatNum','InjType'};
alltransients = exportSessionTransients(data,'sessiontransients_blmin_threshold3SD',analysispath,addvariables);

%% Plot session bin traces with detected transients for each file
% Use plotTransientBins to plot session bins with detected transients for each file.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    allbins = plotTransientBins(data,eachfile,'sigfiltz_normsession_injcropped','sessiontransients_blmin_threshold3SD',maintitle);

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 9, 6]);
    plotfilepath = append(figurepath,'SessionBins_blmin_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end