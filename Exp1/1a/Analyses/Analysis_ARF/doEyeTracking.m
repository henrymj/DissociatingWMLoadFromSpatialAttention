function doEyeTracking(subs)

%quick loop to import edf files to .mat files
subs = [11:12];
root = pwd;
eyeRoot = '/EyeTrack/';
for s = 1:length(subs)
    
    dataRoot = ['/EegData/',char(num2str(subs(s))),'/'];

    sn = subs(s);
    %likely will need to change this for each experiment
    cd([root,eyeRoot])
    %Epoch Data
    dRoot = [root,dataRoot];
    eyeData = Epoch_Data(sn,dRoot);
    save([root,dataRoot,num2str(sn),'_eye_seg.mat'],'eyeData')
    fprintf('Subject Number %d Complete \n',subs(s))
end

cd(root)
