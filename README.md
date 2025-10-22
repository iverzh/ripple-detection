Ripple Detection and Analysis Toolbox

Supplementary code for 
Dickey et al. 2022, Verzhbinsky et al. 2025 - MATLAB scripts for detecting and analyzing hippocampal and neocortical ripple events in neural recordings.
## Overview

This toolbox provides comprehensive methods for:
- Detection of sharp-wave ripple events in broadband LFP recordings
- Cross-correlation analysis between hippocampal and neocortical ripples
- Statistical assessment of ripple coupling and directionality
- Preprocessing of Utah array data

## Requirements

- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Installation

1. Clone the repository:
```bash
git clone https://github.com/iverzh/ripple-detection.git
cd ripple-detection
```

2. Add the code directory to your MATLAB path:
```matlab
addpath(genpath('code/'))
```

## Usage

**Important:** All scripts should be run from the `code/` directory.

### 1. Preprocessing Utah Array Data

```matlab
% Navigate to code directory
cd code/

% Run spike removal and downsampling
RemoveUnitSpikes
```

**Input Requirements:**
- Raw Utah array data files (`.ns5` format)
- Unit spike times structure with fields:
  - `spikeTimes`: Cell array of spike times for each unit
  - `channelID`: Channel identification for each unit
  - `samplingRate`: Sampling frequency of the recording

**Output:**
- Cleaned and downsampled LFP data
- Preprocessed data ready for ripple detection

### 2. Ripple Detection

```matlab
% Run ripple detection on broadband LFP
RunRippleDetection
```

**Input Requirements:**
- Broadband LFP data matrix (channels × time samples)
- Sampling frequency
- Optional: Detection parameters structure with fields:
  - `freqRange`: Frequency range for ripple band (default: [80-250] Hz)
  - `thresholdSD`: Detection threshold in standard deviations
  - `minDuration`: Minimum ripple duration in ms
  - `maxDuration`: Maximum ripple duration in ms

**Output - rippleStats Structure:**
| Field | Description | Format |
|-------|-------------|---------|
| `duration` | Duration of each ripple (ms) | [nRipples × nChannels] |
| `HGz` | Max z-score of high gamma (200+ Hz) amplitude | [nRipples × nChannels] |
| `InterRipPeriod` | Time between adjacent ripples (seconds) | [nRipples-1 × nChannels] |
| `oscFreq` | Peak frequency of each ripple event | [nRipples × nChannels] |
| `rippleAmp` | Maximum analytic amplitude (μV) | [nRipples × nChannels] |
| `window` | Start and end time (samples) | [nRipples × 2 × nChannels] |
| `locs` | Center location of ripples (samples) | [nRipples × nChannels] |
| `locs_sleep` | Ripple locations per data segment | Cell array |
| `density` | Ripple rate (events/minute) | [1 × nChannels] |
| `chanLabel` | Bipolar channel names | Cell array |
| `RejectParams` | Rejection criteria used | Structure |
| `rejectVec` | Rejection condition triggered per event | [nRipples × 1] |

### 3. Cross-Correlation Analysis

#### Hippocampal-Neocortical Analysis
```matlab
% Analyze HC-NC ripple coupling
AnalyzeCrossCorrelogram_HC
```

#### Neocortical-Neocortical Analysis
```matlab
% Analyze NC-NC ripple coupling
AnalyzeCrossCorrelogram_NC
```

**Input Requirements:**
Preprocessed `.mat` file containing `subjPRTH` structure:
- `chanLabels`: Neocortical channel labels (cell array)
- `chanLabelsHC`: Hippocampal channel labels (cell array)
- `HC.eventPRTH`: Cross-correlograms (1ms bins, ±3s window)
  - Format: [nChannelPairs × 6000 bins]
- `HC.nullPRTH`: Null distribution (1000 shuffled iterations)
  - Format: [nChannelPairs × 6000 bins × 1000]
- `HC.nRipples`: Number of ripples per channel pair

**Output Metrics:**
- Percentage of significantly modulated channel pairs
- Percentage showing significant directional coupling
- Percentage where neocortex leads hippocampus
- Statistical significance values (p < 0.05, FDR corrected)
- Cross-correlation peak times and amplitudes

## Advanced Configuration

### Detection Parameters

Customize ripple detection by modifying parameters:

```matlab
params = struct();
params.freqRange = [80 250];        % Ripple frequency band (Hz)
params.thresholdSD = 3;             % Detection threshold (SD above mean)
params.minDuration = 15;            % Minimum duration (ms)
params.maxDuration = 500;           % Maximum duration (ms)
params.minInterRippleInterval = 50; % Minimum time between events (ms)
```

### Null Distribution Settings

Configure statistical testing:

```matlab
nullParams = struct();
nullParams.nShuffles = 1000;        % Number of shuffle iterations
nullParams.shuffleMethod = 'circular'; % 'circular' or 'random'
nullParams.alpha = 0.05;            % Significance level
nullParams.multipleComparison = 'bonferroni'; % Correction method
```

## Troubleshooting

### Common Issues

1. **Memory errors with large datasets:**
   - Process data in segments
   - Reduce sampling rate if appropriate
   - Increase MATLAB heap space

2. **Detection sensitivity:**
   - Adjust `thresholdSD` parameter
   - Verify frequency range matches your ripple characteristics
   - Check signal quality and noise levels

3. **Statistical significance:**
   - Ensure sufficient number of events (>50 recommended)
   - Verify null distribution convergence
   - Consider multiple comparison corrections

## Citation

If you use this code in your research, please cite:

```bibtex
@article{dickey2022,
  title={Cortical Ripples during NREM Sleep and Waking in Humans},
  author={Dickey, Charles W. and Verzhbinsky, Ilya A. and Jiang, Xi and 
          Rosen, Burke Q. and Kajfez, Sophie and Stedelin, Brittany and 
          Shih, Jerry J. and Ben-Haim, Sharona and Raslan, Ahmed M. and 
          Eskandar, Emad N. and Gonzalez-Martinez, Jorge and Cash, Sydney S. and 
          Halgren, Eric},
  journal={Journal of Neuroscience},
  volume={42},
  number={44},
  pages={8378--8389},
  year={2022},
  doi={10.1523/JNEUROSCI.0742-22.2022}
}

@article{verzhbinsky2024co,
  title={Co-occurring ripple oscillations facilitate neuronal interactions between cortical locations in humans},
  author={Verzhbinsky, Ilya A and Rubin, Daniel B and Kajfez, Sophie and Bu, Yiting and Kelemen, Jessica N and Kapitonava, Anastasia and Williams, Ziv M and Hochberg, Leigh R and Cash, Sydney S and Halgren, Eric},
  journal={Proceedings of the National Academy of Sciences},
  volume={121},
  number={1},
  pages={e2312204121},
  year={2024},
  doi={10.1073/pnas.2312204121}
}
```

## License

This project is licensed under the MIT License - see the repository for details.

## Contact

For questions or issues, please open an issue on the [GitHub repository](https://github.com/iverzh/ripple-detection/issues).

## Acknowledgments

This work was supported by NIH grants and tGggOffice of Naval Research. We thank all participants and clinical teams involved in data collection.
