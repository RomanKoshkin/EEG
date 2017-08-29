## Introduction
The Matlab scripts contained in this repository were used to pre-process and analyze the EEG data acquired as part of an ERP study aimed to validate the Efforts Model of Simultaneous Interpreting (Gile, XXXX).The results of the study are presented in a paper that is currently in submission (review).

## Data
The raw and pre-processed EEG and ERP data for all the participants can be found [here](https://cloud.mail.ru/public/3SkP/xwhcS7vXZ). For more information about the dataset structure refer to our paper (Koshkin, Ossadtchi, Shtyrov 2017). 

## Dependencies
* numpy
* Pandas
* Matplotlib

All dependencies can be installed using [pip](https://pip.pypa.io/en/stable/)

## References
[Siraj's video on Dimensionality Reduction](https://www.youtube.com/watch?v=jPmV3j1dAv4&t=281s)

[Detailed answer about PCA on Stats Stackexchange](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues)

[PCA applied to IRIS data set](https://plot.ly/ipython-notebooks/principal-component-analysis/)

[Practical guide to PCA in R and Python](https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/)

## SH
NEW_GRAND_AVERAGE.m 
eeg_event_corr1.m renames event codes (set automatically by the BrainVision Amplifier) to simple numeric values to match the text numbers.
EEG_lag_eval2.m corrects trigger latencies in the EEG dataset to the probes? real onset times based on the audio waveform peaks observed in the BIP channel (#33). Saves the corrected dataset as RAWset2.set.
RAWset3.set ? same as RAWset2 but downsampled from 2000 to 250 Hz, added channel locations, re-referenced to the average mastoid (TP9 TP10), filtered using a zero-phase FIR filter (0.1-30 Hz)

