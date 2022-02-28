Supplementary Files for Dickey et al. 2022

All scrips should be run from the 'code/' directory.

### AnalyzeCrossCorrelogram_HC.m ###
Analyzes the ripple cross correlogram for neocortical ripples centered around hippocampal ripples.
This script computes:
- the percent of channel pairs that are significantly  
- the percent of significantly modulated channel pairs that show signficant sidedness
- the percent of sided channels where the neocortex is leading

The script loads a preprocessed .mat file for each subject that contains a structure (subjPRTH) 
with the following format:

subjPRTH.chanLabels: labels for all neocortical channels
subjPRTH.chanLabelsHC: labels for all hippocamcal channels

subjPRTH.HC.eventPRTH: cross correlograms for each hippocampal-neocortical channel pair in 1 ms 
bins with a window of 6 secs (+- 3 secs)
subjPRTH.HC.nullPRTH: cross correlograms for the null distribution of each channel pair in 1ms. 
Each pair has 1000 iterations of a shuffled null distrubtion to compute significance of the 
real data. 
subjPRTH.HC.nRipples: number of ripples used for each channel pair.

### AnalyzeCrossCorrelogram_NCNC.m ###
Analyzes the ripple cross correlogram for neocortical ripples centered around neocortical ripples in
other locations.
This script computes:
- the percent of channel pairs that are significantly  
- the percent of significantly modulated channel pairs that show signficant sidedness

The format follows the same format as described in the AnalyzeCrossCorrelogram_HC.m section.

### RippleSelection.m ###
This script loads a braodband LFP .mat file (channels x time samples), and runs ripple selection.
The script has two steps:
- Ripple detection, where a set of preliminiary ripple events are detected. This can be turned off 
if these events were already detected/
- Ripple Selection, where the fine tuned rejection criteria are applied to the list of prelim events.

The output rippleStats structure contains the following fields:
- duration: duration of each ripple in ms for each channel.
- HGz: max z score of the high gamma (200+ Hz) analytic amplitude for each ripple.
- InterRipPeriod: period between adjacent ripples within the same channel in seconds.
- oscFreq; the frequency of each ripple event.
- rippleAmp: the maximum analytic amplitude of each ripple in microV.
- window: the start and end time for each ripple in samples.
- locs: locations of the center of each ripple in samples. 
- locs_sleep: locations of the ripples in each data segment (only relevant if multiple
data segments are used);
- density: density / minute of ripple events in each channel.
- chanLabel: the names of each bipolar channel.
- RejectParans: the rejection parameters that were used.
- rejectVec: an array of the rejection condition that was triggered for each preliminary ripple
event.

### PreprocessUtahArray.m ###
This script removes unit spikes and downsamples raw utah array data.
It loads a series of .ns5 files that contain the raw utah array data and an array of unit spike times. 
See data/units.mat for an example of how the unit spike times should be formatted.

