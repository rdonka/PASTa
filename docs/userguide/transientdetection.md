# Transient Detection and Quantification
Transient detection is a critical component of fiber photometry analyses, identifying relevant increases in sensor activation. Previous tools and packages have used a sliding window approach, where all values above an absolute threshold are counted as peaks. Here, we present a novel method of peak detection where each peak is compared to a local baseline and amplitude is calculated and compared to a threshold to determine inclusion. This allows for consistent parameters across the session and reliable detection of individual events despite signal absolute value fluctuation.

## Transient Event Detection
Transients are detected as peaks with a greater amplitude than the specified threshold. To determine amplitude, first a pre-peak baseline must be identified.
PASTa includes three options: window minimum (minimum value within user defined pre-peak window size, e.g. 800ms), window mean (mean of user defined pre-peak window), or local min (absolute local minimum preceding the peak). The amplitude threshold is set by the user, and recommended to be 3SDs. If data are normalized in Z scores, then the criterion is an increase of 3 from baseline. If not, the user inputs the actual value that corresponds to 3SDs in the data stream.

## Transient Event Quantification
Multiple features of transient events can be quantitatively analyzed and compared. Peak detection functions automatically calculate numerous variables for each transient to characterize aspects of both event rise and fall. PASTa includes flexible functions to group transients in the most experimentally relevant manner, such as by time window or experimental condition. 

* __Frequency__: Characterized as peaks per minute. Frequency can be analyzed as whole session frequency, or peaks can be divided into time bins or experimental phases (ITI, during trial, etc).

* __Amplitude__: The height of the event from the pre-peak baseline to the max peak. Note that all events will be at least the value of the set threshold (default 3SD).

* __Rise and Fall Time__: The time difference between the pre-peak baseline and the max peak value. Rise and fall time can also be characterized from half height to peak, allowing analysis of separate rise and fall dynamic shifts.

* __AUC__: Total area under the curve from half height to peak, calculated via the trapezoidal method. Prior to AUC calculation, each transient is linearly transformed so pre-peak baseline is equal to zero.