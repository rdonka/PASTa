# Data Preparation
The first stage in the PASTa Protocol is data preparation, which involves associating fiber photometry data with experimental meta data, loading raw data, and extracting and saving the raw data as MATLAB data structures to faciliate future analysis sessions. For data collected with Tucker Davis Technologies equipment and software _Synapse_, custom functions are included to extract data from raw formats. 

For data collected with other systems, data must be pre-formatted to match a generic csv file format structure and loaded with provided custom functions, or prepared separately by the user. If your data doesn't match the available load options, please feel free to reach out!

## Data Organization 
The PASTa protocol is set up to accomodate any kind of file organization preferred by the user. Organization preference may vary depending on your lab's file storage practices, photometry equipment, or individual projects.

### TDT Synapse Output
Fiber photometry data collected through the software _Synapse_ (Tucker Davis Technologies) is stored in __tanks__ and __blocks__. __Tanks__ are parent folders created by Synapse for each experiment. __Blocks__ are individual folders for each session containing the actual output data. Stored data cannot be accessed directly via the folder, but rather must be extracted via MATLAB.

By default, the tank path is: C:\TDT\Synapse\Tanks. Synapse recognizes experiments and subjects as key categories of information that play a special role in managing data storage and retrieval. Thus, when running the session, it is critical to ensure the correct Experiment profile is selected, and the correct Subject identifier is input for each session. By default, Synapse names data tanks automatically based on experiment name and the start time of the first recording {ExperimentName}-{yymmdd}-{hhmmss}. Blocks of data are named based on subject {SubjectName}-{yymmdd}-{hhmmss} for each recording session and the start time.

Synapse will save a new Tank for every day unless you change the default setting. Click _Menu_ at the top of the bar, then Preferences. Under the _Data Saving_ tab, make sure "New Tank Each Day" is unchecked.

![png](../img/datapreparation_SynapseDataSaving.png)


## Experiment Key Creation
To accomodate a variety of file organization structures, users can first create two csv files containing the information necessary to access raw data, and experimental metadata to match to raw photometry data. MATLAB can access raw data folders stored either locally or in a cloud-based storage app like Box or Dropbox. 

To faciliate analysis, users can create two keys as csv files: a subject key and a file key. In the PASTa protocol, these files will be knit together to pair subject specific information with each individual session of data, preventing the need for manual repeated entry of subject specific information and reduce the time burden of properly maintaining and including experimental metadata factors like subject traits, treatments, and experimental equipment.

First, locate the raw files to be analyzed in your file organization structure. Prepare a folder for extracted raw data to be saved to. This should be separate from your tank folder / raw data storage, and will eventually contain the extracted data sets for each session as MATLAB structures to facilitate efficient data processing in future analysis sessions.

### Subject Key
The subject key should contain information about each subject that is constant and unchanging, such as SubjectID, sex, date of birth, fiber location, sensor, experimental group, and any other user specificed information. 

__For example:__
![png](../img/datapreparation_SubjectKeyExample.png)


### File Key
The file key should contain information about each unique session / file to be analyzed. At a minimum, this must include the SubjectID, folder name, raw folder location, and desired location for the raw data to be exported to.

#### REQUIRED INPUTS

- __SubjectID:__ Unique identifier of the subject. This fieldname should match the first field in the subject key, and inputs for each subject have to match the subject key for subject specific data to be properly associated to fiber photometry data.

- __Folder:__ The name of the folder containing the raw data.

- __RawFolderPath__: The path to the location where the raw data folder is saved for each specific session. All file paths should be specified in the file key WITHOUT the computer user specific portion of the path. This facilitates analysis across multiple devices or cloud storage solutions without manual edits to the file key. Ensure the specified path ends in a forward slash. __For example__, a file saved to a specific device at the path _"C:\Users\rmdon\Box\RawData\"_ should be specified as  _"Box\RawData\"_.

- __ExtractedFolderPath:__ The path to the location where you would like the extracted data structure to be saved. This should be different than the raw data storage location. All file paths should be specified in the file key WITHOUT the computer and user specific portion of the path. This facilitates analysis across multiple devices or cloud storage solutions without manual edits to the file key. Ensure the specified path ends in a forward slash. __For example__, a file saved to a specific device at the path _"C:\Users\rmdon\Box\ExtractedData\"_ should be specified as  _"Box\ExtractedData\"_.

 Any additional fields can be included such as equipment information, recording power, session condition, drug treatments, body weight, and any other variables that are specific to that one session. Note that the only field name that should overlap with a field name in the subject key is _SubjectID_.
    

__For example:__
![png](../img/datapreparation_FileKeyExample.png)

### Making the Experiment Key
The function _loadKeys_ joins the individual subject information to the file key with the data for each session. Additionally, _loadkeys_ appends the unique computer user portion of the file navigation path to the beginning and the Folder name to the end of the raw and extracted folder paths specified in the file key. This creates the full path to the location of each session's data. __For example,__ _"C:\Users\rmdon\Box\RawData\Subject1-240101-121500"_. The created experiment key should be output into a data structure called _experimentkey_.

__Code example:__
![png](../img/datapreparation_loadKeysCode.png)

## Extracting the Data
Prior to beginning analysis, individual session data should be extracted and saved as MATLAB data structures. This makes the process of loading data at the start of each analysis session significantly faster. 

When the raw data is extracted, clipping will be applied by default. This removes the first and last 5 seconds of the session to remove large fluctuations in output signal that occur when the hardware is turned on and off. The number of seconds clipped can be adjusted by overriding the default.

Multiple options are available to extract the data, and should be used depending on the method by which the data was collected, both documented in detail below. Several customized functions are available to extract and format the data for easy processing and analysis. Data collected with TDT can be extracted with custom functions. 


For all other systems, utilize the generic csv format. 


If data were collected with TDT's Synapse, use the function _loadTDTdata_. If data were collected with other systems, use the generic format and load function _loadCSVdata_.

### Set up inputs
Regardless of which extract function you use, the required inputs are the same. 

#### REQUIRED INPUTS
- __rawfolderpaths:__ A string array containing the full paths to the folder locations of raw data to be extracted for each session. This should be formatted as a single column with each full path in a separate row. This is easy to create from the experiment key (see below for an example).

- __extractedfolderpaths:__ A string array containing the full paths to the folder locations where extracted data should be saved for each session (including the individual session name). As with the rawfolderpaths, this should be formatted as a single column with each full path in a separate row and can easily be created from the experiment key (see below for an example).

- __sigstreamnames:__ A cell array containing the strings with the names of all streams to be treated as signal. This allows for flexibility if different photometry rigs have differing naming conventions for the signal stream. Include all signal stream name variations in this cell array. Note that only one stream per file can be treated as signal.

- __baqstreamnames:__ A cell array containing the strings with the names of all streams to be treated as background. This allows for flexibility if different photometry rigs have differing naming conventions for the background stream. Include all background stream name variations in this cell array. Note that only one stream per file can be treated as signal.

- __sigstreamnames:__ A cell array containing the strings with the names of all streams to be treated as signal. This allows for flexibility if different photometry rigs have differing naming conventions for signal stream. Include all signal stream name variations in this cell array. Note that only one stream per file can be treated as signal.

#### OPTIONAL INPUTS
- __clip:__ Number of seconds to clip at the beginning and end of the session. This defaults to 5 seconds.

- __skipexisting:__ This input allows users to toggle if previously extracted raw data files are re-extracted. By default, previously extracted files will be skipped (skipexisting = 1). To override and re-extract all files, set skipexisting = 0.

Extracted raw data files are saved to the extracted folder path as individual MATLAB structures.

## Loading the Data
After raw data is extracted, it can be matched to the experimentkey to associate any subject and session metadata with the photometry data. All sessions will be loaded into one data structure (typically called _data_), to ensure that all sessions are analyzed in the same way throughout following steps of the protocol. Each session of data is a row within the data structure.







This function will use the file key to load each session's extracted data structure and combine them into the structure "data" for analysis.