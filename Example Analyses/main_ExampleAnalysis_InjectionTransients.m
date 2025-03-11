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
    rawdata(eachfile).sessionstart = rawdata(eachfile).injt(1) - (floor(preinjectionlength*60*rawdata(eachfile).fs)); % Find and save start index to rawdata structure field
    sessionend = rawdata(eachfile).injt(2) + (floor(postinjectionlength*60*rawdata(eachfile).fs)); % Find end index and check that it's within the length of the full session
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

[data] = subtractFPdata(data,sigfield,baqfield,fsfield,'artifactremoval',1); % adds sigsub and sigfilt to data frame

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

%% Plot whole session FFT magnitude plots for each file
% Use plotFFTmag to plot all frequency magnitude plots - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    allffts = plotFFTmag(data,eachfile,maintitle,'fs','xmax',100);

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]);
    plotfilepath = append(figurepath,'SessionFFTmag_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Plot whole session FFT power plots for each file
% Use plotFFTpower to plot all frequency magnitude plots - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    allffts = plotFFTpower(data,eachfile,maintitle,'fs');

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 9]);
    plotfilepath = append(figurepath,'SessionFFTpower_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
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

%% Plot whole session normalized streams for each file
% Use plotNormTraces to plot all raw traces - data needs to contain sig, baq, baq_scaled, sigsub, and sigfilt.
streams = {'sigfiltz_normsession', 'sigfiltz_normbaseline'};
streamtitles = {'Normalized to Whole Session', 'Normalized to Baseline'};

for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    normtraces = plotNormTraces(data,eachfile,streams,'fs',maintitle,streamtitles);
    for eachtile = 1:2
        nexttile(eachtile)
        xline(data(eachfile).injt(1),'--','Injection','Color','#C40300','FontSize',8)
    end    

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]);
    plotfilepath = append(figurepath,'SessionNormTraces_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end


%% Find session transients
% Prepare thresholds - since Z scored streams will be analyzed, input threshold as the desired SD.
for eachfile = 1:length(data)
    data(eachfile).threshold3SD = 3;
    data(eachfile).threshold3SDsigfilt = 3*std(data(eachfile).sigfilt_injcropped(1:data(eachfile).injt(1)));
end

% Find session transients based on pre-peak baseline window minimum - reccomended as the first pass choice for transient analysis
[data] = findSessionTransients(data,'blmin','sigfiltz_normbaseline_injcropped','threshold3SD','fs'); % Z scored data (normalized to pre-injection mean)
[data] = findSessionTransients(data,'blmin','sigfilt_injcropped','threshold3SDsigfilt','fs'); % df/f data (threshold calculated as 3*SD of the pre-injection period)

% Find session transients based on pre-peak baseline window mean
[data] = findSessionTransients(data,'blmean','sigfiltz_normbaseline_injcropped','threshold3SD','fs');
[data] = findSessionTransients(data,'blmean','sigfilt_injcropped','threshold3SDsigfilt','fs');

% Find session transients based on pre-peak local minimum (last minumum before the peak in the baseline window)
[data] = findSessionTransients(data,'localmin','sigfiltz_normbaseline_injcropped','threshold3SD','fs');
[data] = findSessionTransients(data,'localmin','sigfilt_injcropped','threshold3SDsigfilt','fs');


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

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 6]);
    plotfilepath = append(figurepath,'SessionBins_blmin_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Plot all transient traces
% Use plotTransientTraces to plot session bins with detected transients for each file.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    tracebins = plotTransientTraces(data,eachfile,maintitle,'fs','sessiontransients_blmin_threshold3SD');

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 6]);
    plotfilepath = append(figurepath,'SessionTransientTraces_blmin_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end


%% Plot transient traces by bin
% Use plotTransientTracesBins to plot session bins with detected transients for each file.
for eachfile = 1:length(data)
    maintitle = append(num2str(data(eachfile).Subject),' - Treatment: ',data(eachfile).InjType); % Create title string for current plot
    tracebins = plotTransientTracesBins(data,eachfile,maintitle,'fs','sessiontransients_blmin_threshold3SD');

    set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 8]);
    plotfilepath = append(figurepath,'SessionBinsTransientTraces_blmin_',num2str(data(eachfile).Subject),'_',data(eachfile).InjType,'.png');
    exportgraphics(gcf,plotfilepath,'Resolution',300)
end

%% Summarize transients by treatment
whichfs = 'fs';
whichsignal = 'sigfiltz_normsession_injcropped';
whichtransients = 'sessiontransients_blmin_threshold3SD';

% Summarize transient means by subject for whole session, pre-injection, and post-injection
subjectsummary_sessiontransients_blmin_threshold3SD = [];
for eachfile = 1:length(data)
    sessionlengthmin = length(data(eachfile).(whichsignal))/(data(eachfile).(whichfs)*60);
    sessionlengthmin_preinjt = data(eachfile).injt(1)/(data(eachfile).(whichfs)*60);
    sessionlengthmin_postinjt = (length(data(eachfile).(whichsignal))-data(eachfile).injt(1))/(data(eachfile).(whichfs)*60);
    
    transients_pre = data(eachfile).(whichtransients).transientquantification.maxloc<data(eachfile).injt(1);
    transients_post = data(eachfile).(whichtransients).transientquantification.maxloc>data(eachfile).injt(1);

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).SubjectID = data(eachfile).SubjectID;
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).InjType = data(eachfile).InjType;

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).freqpermin = height(data(eachfile).(whichtransients).transientquantification) / sessionlengthmin;
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).freqpermin_pre = sum(transients_pre) / sessionlengthmin_preinjt;
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).freqpermin_post = sum(transients_post) / sessionlengthmin_postinjt;

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanamp = mean(data(eachfile).(whichtransients).transientquantification.amp, 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanamp_pre = mean(data(eachfile).(whichtransients).transientquantification.amp(transients_pre), 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanamp_post = mean(data(eachfile).(whichtransients).transientquantification.amp(transients_post), 'omitmissing');

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanrisems = mean(data(eachfile).(whichtransients).transientquantification.risems, 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanrisems_pre = mean(data(eachfile).(whichtransients).transientquantification.risems(transients_pre), 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanrisems_post = mean(data(eachfile).(whichtransients).transientquantification.risems(transients_post), 'omitmissing');

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanfallms = mean(data(eachfile).(whichtransients).transientquantification.fallms, 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanfallms_pre = mean(data(eachfile).(whichtransients).transientquantification.fallms(transients_pre), 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanfallms_post = mean(data(eachfile).(whichtransients).transientquantification.fallms(transients_post), 'omitmissing');

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanAUC = mean(data(eachfile).(whichtransients).transientquantification.AUC, 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanAUC_pre = mean(data(eachfile).(whichtransients).transientquantification.AUC(transients_pre), 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanAUC_post = mean(data(eachfile).(whichtransients).transientquantification.AUC(transients_post), 'omitmissing');

    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanIEIs = mean(data(eachfile).(whichtransients).transientquantification.IEIs, 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanIEIs_pre = mean(data(eachfile).(whichtransients).transientquantification.IEIs(transients_pre), 'omitmissing');
    subjectsummary_sessiontransients_blmin_threshold3SD(eachfile).meanIEIs_post = mean(data(eachfile).(whichtransients).transientquantification.IEIs(transients_post), 'omitmissing');
end


% Summarize transient means by treatment
whichtreatments = {'Saline', 'Morphine'};
overalltransientsummary = [];

for eachtreatment = 1:length(whichtreatments)
    whichfiles = strcmp(char(whichtreatments(eachtreatment)), {subjectsummary_sessiontransients_blmin_threshold3SD.InjType});

    overalltransientsummary(eachtreatment).InjType = char(whichtreatments(eachtreatment)); 

    overalltransientsummary(eachtreatment).freqpermin = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).freqpermin]);
    overalltransientsummary(eachtreatment).freqpermin_pre = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).freqpermin_pre]);
    overalltransientsummary(eachtreatment).freqpermin_post = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).freqpermin_post]);

    overalltransientsummary(eachtreatment).amp = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanamp]);
    overalltransientsummary(eachtreatment).amp_pre = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanamp_pre]);
    overalltransientsummary(eachtreatment).amp_post = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanamp_post]);

    overalltransientsummary(eachtreatment).risems = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanrisems]);
    overalltransientsummary(eachtreatment).risems_pre = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanrisems_pre]);
    overalltransientsummary(eachtreatment).risems_post = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanrisems_post]);

    overalltransientsummary(eachtreatment).fallms = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanfallms]);
    overalltransientsummary(eachtreatment).fallms_pre = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanfallms_pre]);
    overalltransientsummary(eachtreatment).fallms_post = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanfallms_post]);

    overalltransientsummary(eachtreatment).AUC = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanAUC]);
    overalltransientsummary(eachtreatment).AUC_pre = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanAUC_pre]);
    overalltransientsummary(eachtreatment).AUC_post = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanAUC_post]);

    overalltransientsummary(eachtreatment).IEIs = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanIEIs]);
    overalltransientsummary(eachtreatment).IEIs_pre = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanIEIs_pre]);
    overalltransientsummary(eachtreatment).IEIs_post = mean([subjectsummary_sessiontransients_blmin_threshold3SD(whichfiles).meanIEIs_post]);
end

%% Plot overall summary figures
% Bar graph of transient frequency (transients per minute) by treatment condition (Saline vs Morphine)
freqymax = ceil(max([overalltransientsummary.freqpermin, overalltransientsummary.freqpermin_pre, overalltransientsummary.freqpermin_post]) + 1);

close all
frequencyplot = tiledlayout(1,3,'TileSpacing','compact');

nexttile
hold on
freqbar = bar([overalltransientsummary.freqpermin],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
freqbar.CData(1,:) = [0 0 1];  % Blue for Saline
freqbar.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 freqymax]);
ylabel('Transients Per Minute')
title('Transient Frequency (Whole Session)')
hold off

nexttile
hold on
freqbarpre = bar([overalltransientsummary.freqpermin_pre],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
freqbarpre.CData(1,:) = [0 0 1];  % Blue for Saline
freqbarpre.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 freqymax]);
ylabel('Transients Per Minute')
title('Transient Frequency (Pre-Injection)')
hold off

nexttile
hold on
freqbarpost = bar([overalltransientsummary.freqpermin_post],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
freqbarpost.CData(1,:) = [0 0 1];  % Blue for Saline
freqbarpost.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 freqymax]);
ylabel('Transients Per Minute')
title('Transient Frequency (Post-Injection)')
hold off

set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 4]);
plotfilepath = append(figurepath,'OverallTransientBar_Frequency','.png');
exportgraphics(gcf,plotfilepath,'Resolution',300)


% Bar graph of transient amplitude by treatment condition (Saline vs Morphine)
ampymax = ceil(max([overalltransientsummary.amp, overalltransientsummary.amp_pre, overalltransientsummary.amp_post]));

close all
ampplot = tiledlayout(1,3,'TileSpacing','compact');

nexttile
hold on
ampbar = bar([overalltransientsummary.amp],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
ampbar.CData(1,:) = [0 0 1];  % Blue for Saline
ampbar.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 ampymax]);
yline(3);
ylabel('Z Score')
title('Transient Amplitude (Whole Session)')
hold off

nexttile
hold on
ampbarpre = bar([overalltransientsummary.amp_pre],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
ampbarpre.CData(1,:) = [0 0 1];  % Blue for Saline
ampbarpre.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 ampymax]);
yline(3);
ylabel('Z Score')
title('Transient Amplitude (Pre-Injection)')
hold off

nexttile
hold on
ampbarpost = bar([overalltransientsummary.amp_post],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
ampbarpost.CData(1,:) = [0 0 1];  % Blue for Saline
ampbarpost.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 ampymax]);
yline(3);
ylabel('Z Score')
title('Transient Amplitude (Post-Injection)')
hold off

set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 4]);
plotfilepath = append(figurepath,'OverallTransientBar_Amplitude','.png');
exportgraphics(gcf,plotfilepath,'Resolution',300)

% Bar graph of transient risems by treatment condition (Saline vs Morphine)
risemsymax = ceil(max([overalltransientsummary.risems, overalltransientsummary.risems_pre, overalltransientsummary.risems_post]))+50;

close all
risemsplot = tiledlayout(1,3,'TileSpacing','compact');

nexttile
hold on
risemsbar = bar([overalltransientsummary.risems],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
risemsbar.CData(1,:) = [0 0 1];  % Blue for Saline
risemsbar.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 risemsymax]);
ylabel('Rise (ms)')
title('Transient Rise Duration (Whole Session)')
hold off

nexttile
hold on
risemsbarpre = bar([overalltransientsummary.risems_pre],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
risemsbarpre.CData(1,:) = [0 0 1];  % Blue for Saline
risemsbarpre.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 risemsymax]);
ylabel('Rise (ms)')
title('Transient Rise Duration (Pre-Injection)')
hold off

nexttile
hold on
risemsbarpost = bar([overalltransientsummary.risems_post],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
risemsbarpost.CData(1,:) = [0 0 1];  % Blue for Saline
risemsbarpost.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 risemsymax]);
ylabel('Rise (ms)')
title('Transient Rise Duration (Post-Injection)')
hold off

set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 4]);
plotfilepath = append(figurepath,'OverallTransientBar_RiseDuration','.png');
exportgraphics(gcf,plotfilepath,'Resolution',300)

% Bar graph of transient fallms by treatment condition (Saline vs Morphine)
fallmsymax = ceil(max([overalltransientsummary.fallms, overalltransientsummary.fallms_pre, overalltransientsummary.fallms_post]))+50;

close all
fallmsplot = tiledlayout(1,3,'TileSpacing','compact');

nexttile
hold on
fallmsbar = bar([overalltransientsummary.fallms],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
fallmsbar.CData(1,:) = [0 0 1];  % Blue for Saline
fallmsbar.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 fallmsymax]);
ylabel('Fall (ms)')
title('Transient Fall Duration (Whole Session)')
hold off

nexttile
hold on
fallmsbarpre = bar([overalltransientsummary.fallms_pre],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
fallmsbarpre.CData(1,:) = [0 0 1];  % Blue for Saline
fallmsbarpre.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 fallmsymax]);
ylabel('Fall (ms)')
title('Transient Fall Duration (Pre-Injection)')
hold off

nexttile
hold on
fallmsbarpost = bar([overalltransientsummary.fallms_post],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
fallmsbarpost.CData(1,:) = [0 0 1];  % Blue for Saline
fallmsbarpost.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 fallmsymax]);
ylabel('Fall (ms)')
title('Transient Fall Duration (Post-Injection)')
hold off

set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 4]);
plotfilepath = append(figurepath,'OverallTransientBar_FallDuration','.png');
exportgraphics(gcf,plotfilepath,'Resolution',300)

% Bar graph of transient AUC by treatment condition (Saline vs Morphine)
AUCymax = ceil(max([overalltransientsummary.AUC, overalltransientsummary.AUC_pre, overalltransientsummary.AUC_post]))+200;

close all
AUCplot = tiledlayout(1,3,'TileSpacing','compact');

nexttile
hold on
AUCbar = bar([overalltransientsummary.AUC],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
AUCbar.CData(1,:) = [0 0 1];  % Blue for Saline
AUCbar.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 AUCymax]);
ylabel('AUC')
title('Transient AUC (Whole Session)')
hold off

nexttile
hold on
AUCbarpre = bar([overalltransientsummary.AUC_pre],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
AUCbarpre.CData(1,:) = [0 0 1];  % Blue for Saline
AUCbarpre.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 AUCymax]);
ylabel('AUC')
title('Transient AUC (Pre-Injection)')
hold off

nexttile
hold on
AUCbarpost = bar([overalltransientsummary.AUC_post],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
AUCbarpost.CData(1,:) = [0 0 1];  % Blue for Saline
AUCbarpost.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 AUCymax]);
ylabel('AUC')
title('Transient AUC (Post-Injection)')
hold off

set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 4]);
plotfilepath = append(figurepath,'OverallTransientBar_AUC','.png');
exportgraphics(gcf,plotfilepath,'Resolution',300)

% Bar graph of transient IEIs by treatment condition (Saline vs Morphine)
IEIsymax = ceil(max([overalltransientsummary.IEIs, overalltransientsummary.IEIs_pre, overalltransientsummary.IEIs_post]))+2;

close all
IEIsplot = tiledlayout(1,3,'TileSpacing','compact');

nexttile
hold on
IEIsbar = bar([overalltransientsummary.IEIs],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
IEIsbar.CData(1,:) = [0 0 1];  % Blue for Saline
IEIsbar.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 IEIsymax]);
ylabel('Inter Event Interval (s)')
title('Transient IEI (Whole Session)')
hold off

nexttile
hold on
IEIsbarpre = bar([overalltransientsummary.IEIs_pre],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
IEIsbarpre.CData(1,:) = [0 0 1];  % Blue for Saline
IEIsbarpre.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 IEIsymax]);
ylabel('Inter Event Interval (s)')
title('Transient IEI (Pre-Injection)')
hold off

nexttile
hold on
IEIsbarpost = bar([overalltransientsummary.IEIs_post],'FaceColor','flat');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Saline', 'Morphine'});
IEIsbarpost.CData(1,:) = [0 0 1];  % Blue for Saline
IEIsbarpost.CData(2,:) = [1 0 0];  % Red for Morphine
xlim([0.5 2.5]);  % Or [0.7 2.3], etc., as desired
ylim([0 IEIsymax]);
ylabel('Inter Event Interval (s)')
title('Transient IEI (Post-Injection)')
hold off

set(gcf, 'Units', 'inches', 'Position', [0, 0, 12, 4]);
plotfilepath = append(figurepath,'OverallTransientBar_IEIs','.png');
exportgraphics(gcf,plotfilepath,'Resolution',300)