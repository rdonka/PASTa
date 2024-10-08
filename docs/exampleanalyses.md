To help users get started with the PASTa Protocol, here some example analyses are included with raw data to walk users through each step of the process.

# Injection Transients
This example is of whole session transient analysis to determine changes in dopaminergic transients following morphine injection.

## Data Summary
Example data are provided from two subjects. Each subject has two fiber photometry sessions conducted on consecutive days: saline and morphine (10mg/kg). Each session consists of a fifteen minute baseline, i.p. injection marked by epocs sent through Med Associates equipment and stored by Synapse, and a 60 minute post injection recording period.

The example analysis script is located in the repository under the folder Example Analyses. Data are available to download on Box here. The subfolder Injection Transients is prepared for the analysis with folders created for extracted data, analysis, and figure outputs. Raw data blocks collected via Synapse (Tucker Davis Technologies) are nested in the Raw Data folder.

If you have any questions or run into problems accessing the files, please feel free to reach out.

## Data Preparation
The first section of the script sets up paths and analysis key inputs. To enable users to switch computers easily, paths are created without the computer and user specific portion. The computer user specific portion of the path is input to the variable computeruserpath and appended to subsequently needed paths. To access files and functions, paths need to be added to MATLAB via the addpath function. genpath is used within addpath to ensure folders and subfolders at the input path are added. Finally, the names of created Subject and File Keys are added to variables subjectkeyname and filekeyname. These keys contain the session specific and subject specific information needed to load the data and analyze the results. For details on keys, see the user guide section on data preparation.

The keys are loaded using loadKeys to join the subject specific data to the session specific information contained in the file key. Paths to raw data blocks and extracted output locations are added to string arrays, which are input to the function extractFPdata to extract and save raw data as MATLAB structs. Extracted structures are loaded by the function loadKeydata. To trim excess samples before the start of the program and after the end of the post injection period, indices are prepared for each file and the function trimFPdata is used to adjust the streams.

## Signal Processing
To control for motion artifacts and photobleaching, the 405 nm channel (baq) is subtracted from the 465 nm channel (sig) with the function subtractFPdata. The function is called twice to output both the subtracted (sigsub) as well as the subtracted and filtered signal (sigfilt). The scaled 405 nm background channel is also output as _baqsub _and _baqfilt _respectively. After subtraction, sigfilt is normalized. To normalize to the entire session, the function normSession is called, outputting both df/f (sig_normsession_df) and z scored signal streams (sig_normsession_z). To normalize to pre injection baseline, the function normBaseline is called, which uses the pre-injection period to normalize the entire session, also outputting df/f (sig_normbl_df) and z scored (sig_normbl_z) streams. Finally, a for loop is used to create trimmed streams with the injection time period removed. For details on the signal processing functions, see Wiki section 4. Signal Processing.

## Transient Detection
To identify transient events, the findSessionTransients functions are used. These functions detect relevant transients based on an amplitude inclusiong criterion (eg, 3SD amplitude) relative to baseline. Multiple baseline options exist - see function documentation for findSessionTransients_blmin, findSessionTransients_blmean, and findSessionTransients_localmin.

After detection, transients are quantified by frequency, amplitude, half height rise time, half height fall time, and halfheight AUC. After detection and quantification, transients are binned into 5 minute increments. To do so, the number of samples per bin is identified and added to the data structure. The function binSessionTransients identifies the bin of each transient event in the table output by findSessionTransients.

To export the transient events for statistical analysis, the table for each file is appended into the table alltransients which is output as a csv to the location specified by analysispath.

## Plots
First, full session traces are plotted for all relevant streams. In this example, the raw 465 (sig), 405 (baq), subtracted signal (sigsub), subtracted and filtered signal (sigfilt), and baseline normalized signal (sig_normbl_z) are plotted. The function plotTraces is used within a for loop by file to generate plots for each session. plotTraces outputs a tiled plot object (alltraces) with each stream plot as a tile. The tiles are modified to add markers for injection start and end. A main title is also added to the full plot. The plot is saved by subject and injection type into the figure folder specified by the figurepath.