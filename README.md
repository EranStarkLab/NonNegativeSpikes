# NonNegativeSpikes
Code used in Someck et al., 2023, Positive and biphasic extracellular waveforms correspond to return currents and axonal spikes

## Analysis
- waveform_classification.m
  - classification of waveforms and units.
 
## Data
- s_exampledata_mDS2_07.mat
  - Example session from neocortex.
  - structure contains
     mean     mean waveform, matrix of n_samples x n_channels
     sd       SD of waveform, matrix of n_samples x n_channels

 ### To run the code
- Download code and data.
- In MATLAB, follow the instructions at the end of the waveform_classification.m file.