# Process EEG Toolkit
Repository for scripts that can be implemented to pre-process, process and post process EEG datasets.
Current scripts:
> EEGpreprocess_segment.m :  
Segment .edf data and save in .mat / .set / .edf format. (Contributors: ZK)
Dependencies:
- edfread.m (included)
- eeglab (download link redirecting to OneDrive included)
- Sample data can be found in 'sample_data' folder. After downloading and [WIN]extracting eeglab, it is ready to use.

> EEGpreprocess_filternprune.m :  
Remove artifacts (band-pass FIR filter, component based filter) from .edf/.set/.mat files, save in .mat / .set / .edf format. (Contributors: ZK)  
Dependencies:
- eeglab (download link redirecting to OneDrive included)
- (optional) ADJUST: Artifact remover plugin for EEGLab. Currently not impleneted due to source error.
- (optional) MARA: Artifact remover plugin for EEGLab.
- (optional) EEGHEADER_FAUX.mat: Template 'EEG' file for .mat injection into EEGLab. Requires most aspects of the faux and real file to match, so only use .mat loading if necessary {UNTESTED}
- Sample data is created by using EEGpreprocess_segment.m or may be used on the sample data directy via .edf loading.
