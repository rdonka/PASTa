%% EXAMPLE ANALYSIS
% PROJECT SUMMARY: This is an example analysis to determine differences
% in dopamine release to varying frequencies of brain stimulation reward. Recordings are made using GRABDA2H in 
% the nucleus accumbens lateral shell via GRABDA2H. Each session is 30
% minutes with XX trials of intraoral delivery. Within session, all
% intraoral infusions are of the same solution.

%% Set up paths and analysis keys
% Set up user path inputs
computeruserpath =  'C:\Users\rmdon\'; % Computer user unique portion of file path
analysisfolder = 'Box\PASTaDevFiles\Example Analyses\Noncontingent BSR\Analysis\'; % Folder to output analysis csv files to
figurefolder = 'Box\PASTaDevFiles\Example Analyses\Noncontingent BSR\Figures\'; % Folder to output figures to

% Create full paths with computeruserpath appended
analysispath = append(computeruserpath,analysisfolder); 
figurepath = append(computeruserpath,figurefolder);

% Add GitHub Repositories and data folders to MATLAB path
addpath(genpath(append(computeruserpath,'Onedrive\Desktop\GitHub_Repositories\PASTa\'))); % Path for GitHub repository
addpath(genpath(append(computeruserpath,'Box\PASTaDevFiles\Example Analyses\Noncontingent BSR\'))); % Path for analysis files - this is where the keys are saved\\

% Load in experiment key names - Subject Key and File Key
subjectkeyname = 'SubjectKey_ExampleAnalysis_NoncontingentBSR.csv'; % Name of csv file containing subject information; set to '' if not using a Subject Key
filekeyname = 'FileKey_ExampleAnalysis_NoncontingentBSR.csv'; % Name of csv file containing session information and paths

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

extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames,'trim',0,'skipexisting',0); % extract data

%% Load data
% Load previously extracted data blocks and tie to experiment key. 
% Each block is loaded as a row in the data structure.
[rawdata] = loadKeydata(experimentkey); % Load data based on the experiment key into the structure 'rawdata'
 
%% Crop data
% Remove pre and post experimental session samples.
% Use TTL inputs sent by MED for session start and stop - x19__ is start, x20__ is end
% NOTE: Cropping is the only part of the pipeline that will alter the loaded data fields. 
% To ensure you only complete this step once per analysis, it is reccomended to input structure 
% 'rawdata' to the function and output a new structure 'data'.
cropstart = 'x19__'; % name of field with session start index
cropend = 'x20__'; % name of field with session end index
whichstreams = {'sig', 'baq'}; % which streams to crop
whichepocs = {'x19__','x20__', 'x26__','x25__','x27__'}; % which epocs to adjust to maintain relative position - OPTIONAL INPUT

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
    maintitle = append(num2str(data(eachfile).Subject),' - ',data(eachfile).Group); % Create title string for current plot
    alltraces = plotTraces(data,eachfile,maintitle);
    for eachtile = 1:5
        nexttile(eachtile)
        xline(data(eachfile).x26__(1),'--','Block Start','Color','#C40300','FontSize',8)
        for eachblock = 1:length(data(eachfile).x26__)
            xline(data(eachfile).x26__(eachblock),'--','Color','#C40300','FontSize',8)
        end
    end    

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]);
    plotfilepath = append(figurepath,'SessionTraces_',num2str(data(eachfile).Subject),'_',data(eachfile).Experiment,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Plot whole session FFT magnitude plots for each file
% Use plotFFTs to plot all frequency magnitude plots - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - : ',data(eachfile).Group); % Create title string for current plot
    allffts = plotFFTs(data,eachfile,maintitle,'fs');

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]);
    plotfilepath = append(figurepath,'SessionFFTs_',num2str(data(eachfile).Subject),'_',data(eachfile).Experiment,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Normalize data
% To normalize to session mean:
[data] = normSession(data,'sigfilt'); % Outputs whole session z score

% To normalize to a session baseline:
for eachfile = 1:length(data) % prepare indexes for baseline period start and end
    data(eachfile).BLstart = 1;
    data(eachfile).BLend =  data(eachfile).x26__(1); % Start of first block
end
[data] = normBaseline(data,'sigfilt','BLstart','BLend');

%% Plot Normalized Trace
% Use plotNormTrace to plot whole session normalized traces.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - ',data(eachfile).Group, ' - Whole Session Normalization'); % Create title string for current plot
    normtrace = plotNormTrace(data,eachfile,'sigfiltz_normsession','fs',maintitle);

    xline(data(eachfile).x26__(1),'--','Block Start','Color','#C40300','FontSize',8)
    for eachblock = 1:length(data(eachfile).x26__)
        xline(data(eachfile).x26__(eachblock),'--','Color','#C40300','FontSize',8)
    end 

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 3]);
    plotfilepath = append(figurepath,'SessionNormTrace_sigfiltz_normsession_',num2str(data(eachfile).Subject),'_',data(eachfile).Experiment,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

% Use plotNormTrace to plot session baseline normalized traces.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - ',data(eachfile).Group, ' - Session Baseline Normalization'); % Create title string for current plot
    normtrace = plotNormTrace(data,eachfile,'sigfiltz_normbaseline','fs',maintitle);

    xline(data(eachfile).x26__(1),'--','Block Start','Color','#C40300','FontSize',8)
    for eachblock = 1:length(data(eachfile).x26__)
        xline(data(eachfile).x26__(eachblock),'--','Color','#C40300','FontSize',8)
    end 

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 3]);
    plotfilepath = append(figurepath,'SessionNormTrace_sigfiltz_normbaseline_',num2str(data(eachfile).Subject),'_',data(eachfile).Experiment,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end


%% Find session transients
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
end

% Find session transients based on pre-peak baseline window minimum - reccomended as the first pass choice for transient analysis
[data] = findSessionTransients(data,'blmin','sigfiltz_normsession','threshold3SD','fs');

% Bin session transients
[data] = binSessionTransients(data,'sigfiltz_normsession','fs','sessiontransients_blmin_threshold3SD');


% Export transients with added fields for subject and treatment using the EXPORTSESSIONTRANSIENTS function
addvariables = {'Subject','TreatNum','InjType'};
alltransients = exportSessionTransients(data,'sessiontransients_blmin_threshold3SD',analysispath,addvariables);

%% Plot session bin traces with detected transients for each file
% Use plotTransientBins to plot session bins with detected transients for each file.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    allbins = plotTransientBins(data,eachfile,'sigfiltz_normsession_injcropped','sessiontransients_blmin_threshold3SD',maintitle);

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 6]);
    plotfilepath = append(figurepath,'SessionBins_blmin_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Plot transients - overlaid by bin
