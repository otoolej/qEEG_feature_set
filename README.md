[![arXiv](https://img.shields.io/badge/arXiv-1704.05694%20-red)](https://arxiv.org/abs/1704.05694)


NEURAL: A Neonatal EEG Feature Set in Matlab 
============================================

Matlab code to generate a set of quantitative features from multichannel EEG
recordings. Features include amplitude measures, spectral measures, and basic connectivity
measures (across hemisphere's only). Also, for preterm EEG (assuming gestational age < 32
weeks), will generate features from bursts annotations (e.g. maximum inter-burst
interval). Burst annotations require a separate package, also available on
[github](https://github.com/otoolej/burst_detector). 

Releases archived at Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1052810.svg)](https://doi.org/10.5281/zenodo.1052810)


Full details of the methods are in:  

`JM O’Toole and GB Boylan (2017). NEURAL:
quantitative features for newborn EEG using Matlab. ArXiv e-prints, arXiv:1704.05694`  

which is available at https://arxiv.org/abs/1704.05694.

---

[Requirements](#requirements) | [Use](#use) | [Quantitative
features](#quantitative-features) | [Files](#files) | [Test computer
setup](#test-computer-setup) | [Licence](#licence) | [References](#references) |
[Contact](#contact)

## Requirements
Matlab (R2013 or newer,
[Mathworks](http://www.mathworks.co.uk/products/matlab/)) with the signal
processing toolbox and statistics toolbox. Not tested with Octave but should work with minor
tweaking.

## Use 

Set paths in Matlab, or do so using the `load_curdir` function:
```matlab
  >> load_curdir;
```

As an example, generate simulated EEG and calculate relative spectral power, standard
deviation of range-EEG, and brain symmetry index:
```matlab
	% generate EEG-like data (coloured Gaussian noise)
	data_st=gen_test_EEGdata(5*60,64,1);

	% define feature set (or can define in neural_parameters.m):
	feature_set={'spectral_relative_power','rEEG_SD', 'connectivity_BSI'};
	
	% estimate features:
	feat_st=generate_all_features(data_st,[],feature_set);
```

See the `demos/` directory for further examples. All parameters are set the file
`neural_parameters.m`.


## Quantitative features

The feature set contains amplitude, spectral, connectivity, and burst annotation features.
Amplitude features include range-EEG (D. O’ Reilly et al., 2012;
see [references](#references)), a clearly-defined alternative to amplitude-integrated EEG
(aEEG). All features are generated for four different frequency bands (typically 0.5–4,
4–7, 7–13, and 13–30 Hz), with some exceptions. The following table describes the features
in more detail:

| feature name               | description                                                                   | FB  |
|----------------------------|-------------------------------------------------------------------------------|-----|
| spectral\_power            | spectral power: absolute                                                      | yes |
| spectral\_relative\_power  | spectral power: relative (normalised to total spectral power)                 | yes |
| spectral\_flatness         | spectral entropy: Wiener (measure of spectral flatness)                       | yes |
| spectral\_entropy          | spectral entropy: Shannon                                                     | yes |
| spectral\_diff             | difference between consecutive short-time spectral estimates                  | yes |
| spectral\_edge\_frequency  | cut-off frequency (fc): 95% of spectral power contained between 0.5 and fc Hz | no  |
| FD                         | fractal dimension                                                             | yes |
| amplitude\_total\_power    | time-domain signal: total power                                               | yes |
| amplitude\_SD              | time-domain signal: standard deviation                                        | yes |
| amplitude\_skew            | time-domain signal: skewness (absolute value)                                 | yes |
| amplitude\_kurtosis        | time-domain signal: kurtosis                                                  | yes |
| amplitude\_env\_mean       | envelope: mean value                                                          | yes |
| amplitude\_env\_SD         | envelope: standard deviation (SD)                                             | yes |
| connectivity\_BSI          | brain symmetry index (see Van Putten 2007)                                    | yes |
| connectivity\_corr         | correlation (Spearman) between envelopes of hemisphere-paired channels        | yes |
| connectivity\_coh\_mean    | coherence: mean value                                                         | yes |
| connectivity\_coh\_max     | coherence: maximum value                                                      | yes |
| connectivity\_coh\_freqmax | coherence: frequency of maximum value                                         | yes |
| rEEG\_mean                 | range EEG: mean                                                               | yes |
| rEEG\_median               | range EEG: median                                                             | yes |
| rEEG\_lower\_margin        | range EEG: lower margin (5th percentile)                                      | yes |
| rEEG\_upper\_margin        | range EEG: upper margin (95th percentile)                                     | yes |
| rEEG\_width                | range EEG: upper margin - lower margin                                        | yes |
| rEEG\_SD                   | range EEG: standard deviation                                                 | yes |
| rEEG\_CV                   | range EEG: coefficient of variation                                           | yes |
| rEEG\_asymmetry            | range EEG: measure of skew about median                                       | yes |
| IBI\_length\_max           | burst annotation: maximum (95th percentile) inter-burst interval              | no  |
| IBI\_length\_median        | burst annotation: median inter-burst interval                                 | no  |
| IBI\_burst\_prc            | burst annotation: burst percentage                                            | no  |
| IBI\_burst\_number         | burst annotation: number of bursts                                            | no  |

FB: features generated for each frequency band (FB)


## Files
Some Matlab files (.m files) have a description and an example in the header. To read this
header, type `help <filename.m>` in Matlab.  Directory structure is as follows: 
```
├── amplitude_features/                  # amplitude features
├── spectral_features/                   # spectral features
├── connectivity_features/               # hemisphere connectivity features
├── range_EEG/                           # range-EEG (similar to aEEG)
├── IBI_features/                        # features from the burst annotations 
├── preprocessing/                       # loads EEG from EDF files (including artefact removal)
└── utils/                               # misc. functions
```
with some files of interest:
```
├── neural_parameters.m                  # all parameters defined here
├── all_features_list.m                  # complete list of functions (do not edit)
└── generate_all_features.m              # main function: generates feature set on EEG
```



## Test computer setup
- hardware:  Intel Core i7-8700K @ 4.5GHz; 32GB memory.
- operating system: Ubuntu GNU/Linux x86_64 distribution (19.10) with Linux kernel
  5.3.0-24-generic
- software: Matlab (R2019a)


## Licence

```
Copyright (c) 2016, John M. O' Toole, University College Cork
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

  Neither the name of the University College Cork nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```


## References

1. JM O’Toole and GB Boylan (2017). NEURAL: quantitative features for newborn EEG using
Matlab. ArXiv e-prints, arXiv:[1704.05694](https://arxiv.org/abs/1704.05694).

2. D O’Reilly, MA Navakatikyan, M Filip, D Greene, & LJ Van Marter (2012). Peak-to-peak
amplitude in neonatal brain monitoring of premature infants. Clinical Neurophysiology,
123(11):2139–2153.

3. MJAM van Putten (2007). The revised brain symmetry index. Clinical Neurophysiology,
118(11):2362–2367.

4. T Higuchi (1988). Approach to an irregular time series on the basis of the fractal theory,
Physica D: Nonlinear Phenomena, 31:277–283.

5. MJ Katz (1988). Fractals and the analysis of waveforms. Computers in Biology and
Medicine, 18(3):145–156.

6. AV Oppenheim, RW Schafer. Discrete-Time Signal Processing. Prentice-Hall, Englewood
Cliffs, NJ 07458, 1999.

7. JM O’ Toole, GB Boylan, S Vanhatalo, NJ Stevenson (2016). Estimating functional brain
maturity in very and extremely preterm neonates using automated analysis of the
electroencephalogram. Clinical Neurophysiology, 127(8):2910–2918

8. JM O’ Toole, GB Boylan, RO Lloyd, RM Goulding, S Vanhatalo, NJ Stevenson
    (2017). Detecting Bursts in the EEG of Very and Extremely Premature Infants using a
    Multi-Feature Approach. Medical Engineering and
    Physics, vol. 45, pp. 42-50, 2017. 
	[DOI:10.1016/j.medengphy.2017.04.003](https://dx.doi.org/10.1016/j.medengphy.2017.04.003)
	
9. JM O'Toole and GB Geraldine (2019). Quantitative Preterm EEG Analysis: The Need for
   Caution in Using Modern Data Science Techniques. Frontiers in Pediatrics 7, 174
   [DOI:10.3389/fped.2019.00174](https://doi.org/10.3389/fped.2019.00174)
	

## Contact

John M. O'Toole

Neonatal Brain Research Group,  
[INFANT Research Centre](https://www.infantcentre.ie/),  
Department of Paediatrics and Child Health,  
Room 2.19 UCC Academic Paediatric Unit, Cork University Hospital,  
University College Cork,  
Ireland

- email: jotoole AT ucc _dot ie 

