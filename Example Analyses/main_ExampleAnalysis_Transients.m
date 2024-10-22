%% EXAMPLE ANALYSIS
% PROJECT SUMMARY: This is an example analysis to determine changes in 
% GCaMP6f transients in VTA dopamine neurons after within session saline or 
% morphine injection. Each subject has two recording sessions: a saline 
% control session and a morphine administration session. Each session
% recording consists of a 15 minute pre-injection baseline, injection, and 
% a 60 minute post injection period..

%% Set up paths and analysis keys
% Set up user path inputs
computeruserpath =  'C:\Users\rmdon\'; % Computer user unique portion of file path
analysisfolder = 'Box\PASTaExampleFiles\Analysis\'; % Folder to output analysis csv files to
figurefolder = 'Box\PASTaExampleFiles\Figures\'; % Folder to output figures to

% Create full paths with computeruserpath appended
analysispath = append(computeruserpath,analysisfolder); 
figurepath = append(computeruserpath,figurefolder);

% Add GitHub Repositories and data folders to MATLAB path
addpath(genpath(append(computeruserpath,'Onedrive\Desktop\GitHub_Repositories\PASTa\'))); % Path for GitHub repository
addpath(genpath(append(computeruserpath,'Box\PASTaExampleFiles\'))); % Path for analysis files - this is where the keys are saved\\

% Load in experiment key names - Subject Key and File Key
subjectkeyname = 'Subject Key - Example Analysis - Morphine Transients.csv'; % Name of csv file containing subject information; set to '' if not using a Subject Key
filekeyname = 'File Key - Example Analysis - Morphine Transients.csv'; % Name of csv file containing session information and paths

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

% Extract data
extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames); % extract data

%% Load data
% Load previously extracted data blocks and tie to experiment key. 
% Each block is loaded as a row in the data structure.
[rawdata] = loadKeydata(experimentkey); % Load data based on the experiment key into the structure 'rawdata'

%% Trim data
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

% TRIM: Trim pre session and post session samples
% Trimming is the only part of the pipeline that will alter the loaded
% data fields. To ensure you only complete this step once per analysis, it
% is reccomended to input structure 'rawdata' to the function and output a
% new structure 'data'.
trimstart = 'sessionstart'; % name of field with session start index
trimend = 'sessionend'; % name of field with session end index
whichstreams = {'sig', 'baq'}; % which streams to trim
whichepocs = {'injt','sess'}; % which epocs to adjust to maintain relative position

[data] = trimFPdata(rawdata,trimstart,trimend, whichstreams,whichepocs); % Output trimmed data into new structure called data

%% Process data
% Subtract and filter data with default settings
sigfield = 'sig';
baqfield = 'baq';
fs = data(1).fs;

[data] = subtractFPdata(data,sigfield,baqfield,fs); % adds sigsub and sigfilt to data frame

% END OF FINALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize data
% To normalize to session mean:
[data] = normSession(data,'sigfilt'); % Outputs whole session z score

% To normalize to a session baseline:
for eachfile = 1:length(data) % prepare indexes for baseline period start and end
    data(eachfile).BLstart = 1;
    data(eachfile).BLend =  data(eachfile).injt(1);
end
[data] = normBaseline(data,'sigfilt','BLstart','BLend','both'); % 'both' will output df/f and z score

%% Remove injection time window from normalized signal
% This loop removes samples between the start and end of the injection. For
% trimmed data, injt(1) will be used as the injection time point.
normstreams = {'sig_normsession_df', 'sig_normsession_z','sig_normbl_df','sig_normbl_z'};

for eachfile = 1:length(data)
    for eachstream = 1:length(normstreams)
        currstream = char(normstreams(eachstream)); % Pull stream name from cell array into a character object

        allindices = (1:length(data(eachfile).(currstream))); % Temporary list of all indices in the stream
        includeindices = (allindices < data(eachfile).injt(1) | allindices > data(eachfile).injt(2)); % Find indices before injection start and after injection end
        data(eachfile).(append(currstream,'_trimmed')) = data(eachfile).(currstream)(includeindices); % Add new fields ending in '_trimmed' to data structure
    end
end

%% Find session transients
% Find all transient events within the session using default parameters
for eachfile = 1:length(data)
    data(eachfile).transthreshold = 3;
end
[data] = findSessionTransients(data,'sig_normbl_z_trimmed');

[data] = findSessionTransients_DEV(data,'sig_normbl_z_trimmed','transthreshold','fs');

%% Bin session transients
% Set up samples per bin
minsperbin = 5; % Length of each bin in minutes
for eachfile = 1:length(data)
    data(eachfile).binsamples = floor(data(eachfile).fs*60*minsperbin); % Bin length in samples using sampling rate fs
end

whichbinsamples = 'binsamples'; % Name of field containing the calculated number of samples per bin
whichstream = 'sig_normbl_z_trimmed'; % Which stream was used for transient analysis
whichtransients = 'sessiontransients'; % Name of field containing the output table of transients from findSessionTransients
whichpklocs = 'pklocs'; % Name of field containing the transient peak locations in the output table from findSessionTransients

[data] = binSessionTransients(data,whichbinsamples,whichstream,whichtransients,whichpklocs);

% Export transients with added fields for subject and treatment
alltransients = table; % Create table of all transients for all files
for eachfile = 1:length(data)
    eachtransients = data(eachfile).sessiontransients;
    eachtransients.Subject(1:size(eachtransients,1)) = {data(eachfile).Subject};
    eachtransients.TreatNum(1:size(eachtransients,1)) = data(eachfile).TreatNum;
    eachtransients.InjType(1:size(eachtransients,1)) = {data(eachfile).InjType};
    alltransients = vertcat(alltransients,eachtransients);
end
writetable(alltransients,append(analysispath,'SalinevsMorphineInjection_TransientsbyBin.csv'));

%% Plot session traces
% Use plotTraces to plot all raw traces - data needs to contain sig, baq,
% sigsub, sigfilt, and normsig_z
whichstreams = {'sig','baq','sigsub','sigfilt','sig_normbl_z'};
whichtitles = {'Raw 465','Raw 405','Subtracted Signal','Subtracted and Filtered Signal','Normalized Signal - Baseline'};
whichylabels = {'f','f','df/f','df/f','z score'};
whichcolors = {'#0E53FF','#600085','#6B39F0','#0DA1A8','#2AB249'};

for eachfile = 1:length(data)
    hold on
    alltraces = plotTraces(data,eachfile,whichstreams,whichtitles,whichylabels,whichcolors)
    for eachtile = 1:length(alltraces)
        xline(alltraces(eachtile),data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',8)
        xline(alltraces(eachtile),data(eachfile).injt(2),'--','Color','#C40300','FontSize',8)
    end
    sgtitle(append('Subject ',num2str(data(eachfile).Subject),' - Treatment ', num2str(data(eachfile).TreatNum)))
    hold off

    finaltraceplot = gcf;
    exportgraphics(finaltraceplot,append(figurepath,'SessionSignalTraces_',data(eachfile).Subject,'Treatment_',data(eachfile).InjType, ...
        '.png'),'Resolution',300)
end

%% Plot session bins with transients
% Use plotTraceBins to plot the session by bin with transients marked. 
% plotTraceBins plots session bins for one file. Use it in a for loop to
% plot and save all of your files.

% Set up inputs
whichbinsamples = 'binsamples'; % Name of field containing the calculated number of samples per bin
whichstream = 'sig_normbl_z_trimmed'; % Which stream was used for transient analysis
whichtransients = 'sessiontransients'; % Name of field containing the output table of transients from findSessionTransients
whichpklocs = 'pklocs'; % Name of field containing the transient peak locations in the output table from findSessionTransients
whichylabel = 'Z Score'; % String with label for the y axis
whichpkcolors = {'#C60069'}; % Color to plot transients - to change color by bin, this can be a cell array the length of the number of bins you have

% Set up color arrays for streams - set by bin to alter color for pre and post injection by treatment. 
% Here we can make two cell arrays and then select which one to apply for each file in the for loop.
streamcolors_saline = {'#0072BD','#0072BD','#0072BD',... % Pre injection (blue)
                        '#0072BD','#0072BD','#0072BD','#0072BD','#0072BD','#0072BD','#0072BD', '#0072BD','#0072BD','#0072BD','#0072BD','#0072BD'}; % Post injection (blue)
streamcolors_morphine = {'#0072BD','#0072BD','#0072BD',... % Pre injection (blue)
                          '#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18','#FFDC18'}; % Post injection (yellow)

% Make bin trace plots
for eachfile = 1:length(data)
    % Set plot stream colors by treatment (blue for saline, yellow for morphine)
    if data(eachfile).TreatNum == 1 % Set plot stream colors by treatment (blue for saline, yellow for morphine)
        whichstreamcolors = streamcolors_saline;
    elseif data(eachfile).TreatNum == 2
        whichstreamcolors = streamcolors_morphine;
    end

    % Plot bins
    hold on
    bintraces = plotTraceBins(data,eachfile,whichbinsamples,whichstream,whichtransients,whichpklocs,whichylabel,whichstreamcolors,whichpkcolors)
    
    % Modify x axis for each subplot tile
    for eachtile = 1:length(bintraces)
        fsmin = floor(data(eachfile).fs*60);
        xticks(bintraces(eachtile),[1 fsmin*1 fsmin*2 fsmin*3 fsmin*4 fsmin*5])
        xticklabels(bintraces(eachtile),{'0','1', '2', '3', '4', '5'})
        xlabel(bintraces(eachtile),'Minute')
    end
    
    % Add an overall title to the figure
    sgtitle(append('Subject ',data(eachfile).Subject,' - ', data(eachfile).InjType))
    hold off
    
    % Save figure
    finalbinplot = gcf;
    exportgraphics(finalbinplot,append(figurepath,'BinTraces_',data(eachfile).Subject,'Treatment_',data(eachfile).InjType, ...
        '.png'),'Resolution',300)
end
