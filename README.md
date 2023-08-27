# NonNegativeSpikes
Code used in Someck et al., 2023, Positive and biphasic extracellular waveforms correspond to return currents and axonal spikes

### Analysis
- waveform_categorization.m
  - Categorization of waveform and unit classification.
 
### Data
- s_exampledata_mDS2_07.mat
  - Example session from neocortex.
  - Structure contains 137 units, for each:
    - nspk:     number of spikes
    - mean:     mean waveform, matrix of n_samples x n_channels
    - sd:       SD of waveform, matrix of n_samples x n_channels
- s_full_dataset.mat
  - The full dataset used in Someck et al., 2023.
  - Structure contains 9160 units, for each:
     - filename: session identifier
     - shankclu: unit identifier [shank number, unit number]
     - nspk:     number of spikes
     - mean:     mean waveform, matrix of n_samples x n_channels
     - sd:       SD of waveform, matrix of n_samples x n_channels
     - region:   neocortex/CA1 0/1
  - Data available at https://doi.org/10.5281/zenodo.8286842
   
 ## To run the code
- Download code and data.
- In MATLAB, follow the instructions at the end of the waveform_categorization.m file.
