function runERPs
% Wrapper function for doERPs
dbstop if error

subs = [38];

root = pwd;
erpRoot =  '/EegPipeline/';

for i = 1:length(subs)
    tic
    % Print subject number
    sn = subs(i);
    fprintf('Subject:\t%d\n',sn)
    fName = [root, '/EegData/', num2str(subs(i)),'/',num2str(subs(i)),'_EEG.mat'];

    eRoot = [root,erpRoot];
    cd(eRoot)
    erp = doERPs(subs(i));

    %we now just save all variables we make as the EEG file
    fprintf('saving file.. \n')
    save(fName,'erp','-v7.3')
    toc

    %%%% don't save unsegmented
    erp.data = [];
    erp.eventTimes = [];
    erp.event = [];
    erp.eventCodes = [];
end

cd(root)
