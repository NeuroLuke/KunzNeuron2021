# KunzNeuron2021
Analysis code for Kunz et al., Neuron, 2021: "A neural code for egocentric spatial maps in the human medial temporal lobe".

# Overview
The repository contains the folders "Figures", "OpenField", and "TreasureHunt".
All code is written in MATLAB (MathWorks, Natick, MA, USA).

# Figures
This folder contains the data and a script to recreate the figures from the paper.
Absolute paths have to be adjusted in order to successfully run this script.

# OpenField
This folder contains the code used to analyze the data from the spatial reference memory task ("OpenField").

The subfolders contain:
- BasicsAboutUnits_300419: Code to calculate quality metrics for the single-neuron data from the OpenField task.
- Beh_210318: Code to process and analyze the behavioral data from the OpenField task.
- CuePresentation_20200803: Code to analyze the activity of object cells, egocentric bearing cells, and object by egocentric bearing cells during spatial memory recall.
- ELRPDCellAnalysis_20200623: Code to analyze egocentric bearing cells ("ELRPD cells").
- ELRPDCellDistanceModulation_20200729: Code to examine distance tuning in egocentric bearing cells.
- ELRPDCellGoalTuning_20201122: Code to examine goal tuning in egocentric bearing cells.
- ELRPDCellTemporalStability_20200727: Code to examine the temporal stability of the tuning of egocentric bearing cells.
- Functions: Functions that support the code in the other subfolders.
- MemoryCellAnalysis_20200807: Code to examine memory cells and their overlap with egocentric bearing cells, direction cells, and place-like cells.
- ObjectCellAnalysis_010120: Code to test for object cells.
- PlaceDirAnalysis_030620: Code to analyze direction cells and place-like cells.
- ReferencePointAnalysis_20200806: Code to examine the reference points.
- SpikeExtraction_NewWC_ManOpt_241219: Code to perform spike detection and sorting.

# TreasureHunt
This folder contains the code used to analyze the data from the hybrid spatial navigation-episodic memory task ("TreasureHunt").

The subfolders contain:
- BasicsAboutUnits_20200902: Code to calculate quality metrics for the single-neuron data from the TreasureHunt task.
- Beh_210420: Code to process the behavioral data from the TreasureHunt task.
- BehAnalysis_280520: Code to analyze the behavioral data from the TreasureHunt task.
- ELRPDCellAnalysis_220420: Code to analyze egocentric bearing cells ("ELRPD cells").
- Functions: Functions that support the code in the other subfolders.
- LocRecall_20200617: Code to analyze the cellular activity during object-cued location recall.
- ObjRecall_20200622: Code to analyze the cellular activity during location-cued object recall.
- PlaceDirAnalysis_220520: Code to analyze direction cells and place-like cells.
- SpikeExtraction_NewWC_ManOpt_20200621: Code to perform spike detection and sorting.
- TowerTransport_20200618: Code to examine the tuning of egocentric bearing cells during passive movement.

# Contact
drlukaskunz@gmail.com

