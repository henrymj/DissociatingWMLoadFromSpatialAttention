function erp = doERPs(sub)
% Preprocesses EEG data.
% All settings to be changed are set in the preprocessSettings function.

clear erp

cur = ['EegPipeline/'];
p.root = pwd;
p.root = p.root(1:end-length(cur));

% Read in data and assign parameters
data_path = ['/Users/Gisella/Documents/MATLAB/Projects/Perceptual Grouping/Perceptual Grouping V 6.0 Sequential/Data/'];
data_file = sprintf('s%02d_PGS.vhdr',sub);
erp = load_BrainProducts_file(data_path,data_file);
erp.data_file = data_file;

% load preprocessing settings
erp = preprocessSettings(erp);

% Re-reference and/or drop unwanted electrodes
erp = Rereference(erp);

% Remove spaces and letters from the BrainProducts codes
erp = ChangeCondCodes(erp); 

% Segment data
erp.trial = segment(erp,erp.preTime,erp.postTime,erp.data); % for full trial
erp.arfDat = segment(erp,erp.arfPreTime,erp.arfPostTime,erp.data); % for artifacts

% Artifact Rejection: mark bad data
erp = artReject(erp);

% Baseline correction
% erp.trial.baselined = doBaseline(erp.trial.data,erp.trial.times,erp.baseStart,erp.baseEnd);  % for full trial
erp.arfDat.baselined = doBaseline(erp.arfDat.data,erp.arfDat.times,erp.baseStart,erp.baseEnd); % for artifacts
