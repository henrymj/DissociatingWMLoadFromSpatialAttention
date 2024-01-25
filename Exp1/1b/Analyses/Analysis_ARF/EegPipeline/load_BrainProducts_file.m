function erp = load_BrainProducts_file(data_path,filename)

fprintf('Importing Data\n');

%% Read the header file
CONF = readbvconf(data_path, filename); %channel info etc.
%% load in the EEG file
[EEG, com]=pop_loadbv(data_path, filename);  % requires EEG lab.
%% Rename things so we have them to use later!!
EEGdata = EEG.data;
erp.data = EEGdata; % changed from eeg to erp to match josh's naming schema
erp.srate = EEG.srate;
erp.nChans = size(EEGdata,1); % Number of Channels
erp.nArfChans = size(EEGdata,1) - 1; % don't count the stim trak as a channel to artifact reject ! 
erp.pnts = size(EEGdata,2); % number of overall sample data points
erp.rateAcq = 1000/EEG.srate; % 2 ms= Rate of Data Acquisition
erp.event = struct( 'type', { EEG.event.type }, 'latency', {EEG.event.latency});
erp.eventTimes = round(cell2mat({erp.event.latency})); % Event Times
erp.eventCodes = nan(1,length(erp.eventTimes)); % Event Codes (S1 S2 R1 R2 etc)  % JJF: This is being preallocated, where does it get used.

erp.headerInfo = CONF.comment;

for c = 1:erp.nChans % loop through channels and get labels and coordintates % JJF: I added this
    
    % grab channel labels
    label = strsplit(CONF.channelinfos{c},',');
    erp.chanLabels{c} = label{1};
       
    % grab channel coordinates
    coordinates = strsplit(CONF.coordinates{c},',');
    erp.chanCoordinates{c} = [str2double(coordinates{1}) str2double(coordinates{2}) str2double(coordinates{3})];

end                            
