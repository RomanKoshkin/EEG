## Introduction

This repository contains MATLAB scripts and links to the data used as part of our EEG/ERP study (Koshkin, Shtyrov & Ossadtchi, 2017) to explore how WM load affects attention to task-irrelevant probe stimuli during simultaneous interpreting. The draft of the paper presenting the key findings is available at [bioRxiv.org](https://cloud.mail.ru/public/3Lon/PVvjRVfY4). Note that this repo contains only MATLAB scripts we used mainly for EEG dataset processing. For other scripts and datasets, refer to our other repositories: [R scripts](https://github.com/RomanKoshkin/R), [VBA](https://github.com/RomanKoshkin/VBA) and [NLP-Python](https://github.com/RomanKoshkin/NLP-Python).

## Scripts and Functions Contained in this Repo

### NEW_GRAND_AVERAGE.m 

This MATLAB script plots ERPs specific to three distinct levels of WM load: low, medium and high. Unlike the limited functionality of the script SUPER.m, this one allow one to easily manipulate different parameters. For example, set subject- and direction-specific boundaries. It also plots effect size heatmaps and interaction plots, as well as boxplots of subject- and direction-specific WM load distributions. Finally, it performs some exploratory statistical tests. The actual tests, however, were performed in R (the scripts are available [here](https://github.com/RomanKoshkin/R)).

### SUPER.m

This MATLAB script is a wrapper for many of the EEG pre-processing routines. Specifically, it
-	gets files from the directories corresponding to subject codes (e.g. KOS_EEG in the project home folder)
-	performs artifact removal using the [Artifact Subspace Reconstruction]( https://www.google.ru/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0ahUKEwjcopyLpP3VAhUJJ5oKHSRzDFkQFggvMAE&url=https%3A%2F%2Fsccn.ucsd.edu%2Feeglab%2Fplugins%2FASR.pdf&usg=AFQjCNEWMNv8JmmqhpschRTKlwo2Lffasg) algorithm;
-	calls the interp_chan.m function as necessary to interpolate missing channels. In some subjects the temporal channels (e.g. FT9 FT10 T7 T8) were very noisy and had to be removed and re-interpolated.
-	changes event codes to match WM load estimates by calling the EVS_corr_SUPER.m function. In turn, EVS_corr_SUPER.m also linearly interpolates WM load values at probe onset latencies. Finally, uncommenting line 79 in EVS_corr_SUPER.m, you will call the plotinterp.m function that plots curves for interpolated and noninterpolated WM load against the time-coded source text transcript and its corresponding time-coded translation.
-	epochs the continuous EEG data using the EPOCH_SUPER.m function.
-	plots ERPs by subject, direction of interpretation and WM load. Unlike the NEW_GRAND_AVERAGE.m function (see below), it does not compute subject- and direction-specific cutoffs. The cutoff values for every WM load estimator are preset in the params.mat file;
-	optionally runs a genetic algorithm (implemented from scratch in the function GA_super.m) to reject “bad” epochs. This algorithm is just a proof of concept and was NOT used for actual data cleaning.
-	calls grand_avg_SUPER.m to plot grand average ERPs corresponding to high and low WM load defined based on cutoffs set in the params.m file.

### EEG_event_corr1.m

Renames event codes (set automatically by the BrainVision amplifier) to simple numeric values to match the text numbers.

### ExtractDataFromExcelPY.m 

This script extracts relevant WM load data and corresponding time codes from the Excel spreadsheets using the MATLAB-Python interface.

### EEG_lag_eval2.m 

Corrects trigger latencies in the EEG dataset to the real onset times of the probes based on the audio waveform peaks observed in the BIP channel (#33). Saves the corrected dataset as RAWset2.set.

### RAWset3.set 

Same dataset as RAWset2, but downsampled from 2000 to 250 Hz, has channel locations, is re-referenced to the average mastoid (TP9 TP10), and filtered using a zero-phase FIR filter (0.1-30 Hz).

### Filt_SUPER.m

This script batch-filters the datasets specified using a zero-phase FIR filter, runs ICA (binica) on the RAWset3_[subject].set, and saves the ICA matrix in the current subject?s directory as ICA.mat, so that there is no need to re-compute the decomposition each time you run the entire pre-processing pipeline.

### myplotX.m

This function takes epoch numbers as the only parameter. It finds these specified epochs in the current EEGlab dataset, averages them and plots the average ERP waveform. The function is convenient for quickly exploring the different combination of epochs. Speeds up time compared to using EEGlab's GUI.

### ICAimport.m

The script loads previously computed ICA information and pastes it into the current dataset. It also opens the dialog for reviewing the independent components, their topographies and marking them for rejection.

### Rand_test.m

This script performs the randomization test of the mean voltage difference between the high and low WM load conditions in sliding windows, on grand average ERPs.

### WMloadDynamics.m

This script computes how average WM load changes over time and creates barplots. Insert a picture

### latency.m

This script performs a randomization test of equal latencies in the windows of interest. The null hypothesis is that the latencies are the same

### GA_SUPER.m

Proof-of-concept algorithm for rejecting “bad” epochs. ![Alt text](https://github.com/RomanKoshkin/EEG/blob/master/GA.png)

### max_projection_fork1.m

This script attempts to find the optimal spatial filter to maximize the difference between low and high WM load in the window of interest.

All the above scripts ran trouble-free under EEGlab v13.5.4b. We cannot guarantee that they would work with other versions of EEGlab.


## Data

Download the datasets [here](https://cloud.mail.ru/public/4cRw/2z7aTxQTx)
