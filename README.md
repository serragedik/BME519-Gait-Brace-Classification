# BME519-Gait-Brace-Classification
This repository contains MATLAB scripts for preprocessing Rokoko gait data, extracting gait-cycle features, and preparing machine-learning tables for brace condition classification.

Files

preprocessing_and_segmentation.m
Processes raw Rokoko CSV files. It performs Butterworth filtering, gait-cycle segmentation using vertical foot position, cycle normalization, quality control, and saves segmented cycles as MAT files.

features.m
Loads the segmented MAT files and extracts time-domain and harmonic-domain gait features. It also prepares the tables used for MATLAB Classification Learner.

How to Run

Place the raw Rokoko CSV files inside subject folders.

Example folder structure:

Project_Folder
Subject_01
trial1.csv
trial2.csv
Subject_02
trial1.csv
trial2.csv

Open MATLAB and run:

preprocessing_and_segmentation

If the folder path inside the script is not valid, MATLAB will ask you to select the root folder containing the subject folders.

After preprocessing, segmented MAT files will be saved under:

_processed_segments/segments_mat/

Run:

features

If the segment folder path inside the script is not valid, MATLAB will ask you to select the folder containing the segmented MAT files.

The feature extraction script creates the expanded feature table and train-test tables for Classification Learner.

Main output files:

cycle_features_for_ml_expanded.csv
cycle_features_for_ml_expanded.mat
E0_H3_expanded_trial_level_TRAIN.csv
E0_H3_expanded_trial_level_TEST.csv
E0_H3_expanded_subject_level_TRAIN.csv
E0_H3_expanded_subject_level_TEST.csv

Machine Learning

Machine-learning models were trained using MATLAB Classification Learner.

Response variable:

braceLabel

Main predictor variables:

cadence, maxFlexion, strideDuration, romFlexion, meanFlexion, minFlexion, stdFlexion, rmsFlexion, harmonic1Amp, harmonic2Amp, harmonic3Amp, spectralEnergy, spectralEntropy, logLowHighHarmonicRatio

Metadata columns such as subjectId, trialId, sourceFile, braceLevel, inclineLevel, and speedLevel were not used as predictors.

Experiments

The generated feature tables were used for:

All-data random validation
Fixed speed-incline random validation
Fixed speed-incline trial-level split
Fixed speed-incline subject-level split

For trial-level and subject-level experiments, models were selected based on validation performance on the training set and then evaluated on held-out test data.

Notes

Raw CSV files and segmented MAT files are not included due to file size limitations.

Model performance was evaluated using accuracy, confusion matrices, ROC curves, precision, recall, and F1-score.
