Reading .nex files
------------------

Use readNexFile function to read the contents of *.nex file in Matlab.
For example:

nexFileData = readNexFile(filePath);
nexFileData1 = readNexFile(); % if file path is not specified, Open File dialog will be used

nexFileData is a structure that contains all the data from .nex file.
See comments in the readNexFile.m file for more information on the contents of nexFileData structure.

For .nex5 files, use readNex5File function:
nex5FileData = readNex5File(filePath);


Writing .nex and .nex5 files
------------------

Use writeNexFile function to write the contents of nexFileData structure to .nex file.
Use write5NexFile function to write the contents of nexFileData structure to .nex file.

You can read .nex file and then save nexFileData in .nex5 file and vice versa.

For example (see file exampleSaveDataInNexAndNex5File.m):

% start new nex file data
nexFile = nexCreateFileData(40000);

% add continuous variable
% digitizing frequency 1000 Hz
Fs = 1000;
% time interval from 1 to 5
t= 1:1/Fs:5;
% sin with frequency 2 Hz
contValues = sin(2*pi*t*2);
% specify start time (t(1)), digitizing frequency (Fs), data (x2) and name
nexFile = nexAddContinuous(nexFile, t(1), Fs, contValues, 'sin2Hz');

% add continuous variable with 2 fragments: 
firstFagmentValues = cos(2*pi*t*2);
secondFagmentValues = sin(2*pi*t*2);

% start of first fragment at 1.5 seconds
% start of second fragment at 10 seconds
fragmentTimestamps = [1.5; 10];

% index of the first element in the first fragment is 1
% index of the first element in the second fragment is 
% 1 + length(firstFagmentValues)
fragmentIndexes = [ 1; 1+length(firstFagmentValues)];
% merge data for all fragments into one vector
allValues = [firstFagmentValues, secondFagmentValues];
nexFile = nexAddContinuousWithMultipleFragments(nexFile, fragmentTimestamps, fragmentIndexes, Fs, allValues, 'twoFragments');

% add neuron spike train
neuronTs = [0.5 0.9 2.1 2.3 2.5]';
nexFile = nexAddNeuron(nexFile, neuronTs, 'neuron1');

% add event spike train
eventTs = [10 20 30 40]';
nexFile = nexAddEvent(nexFile, eventTs, 'event1');

% add interval variable
intStarts = [5 10];
intEnds = [6 12];
nexFile = nexAddInterval(nexFile, intStarts, intEnds, 'interval1');

% add  waveforms
% waveform timestamps
waveTs = [1 2]';
% 2 waveforms (columns), 5 data points each
waves = [-10 0 10 20 30; -15 0 15 25 15]';
nexFile = nexAddWaveform(nexFile, 40000, waveTs, waves, 'wave1');

% save nex file in user's documents
userDir= getenv('USERPROFILE');

nexFilePath = strcat(userDir, '\Documents\test1.nex');
writeNexFile(nexFile, nexFilePath );

% save .nex5 file also
nex5FilePath = strcat(userDir, '\Documents\test1.nex5');
writeNex5File(nexFile, nex5FilePath);

% save .nex5 file with continuous data written as floats
nex5FilePathFloats = strcat(userDir, '\Documents\test1WithFloats.nex5');
writeNex5File(nexFile, nex5FilePathFloats, 1);