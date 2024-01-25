pList = [1,2,4:9,11:12];

root = pwd;
behRoot =  'Data/';

for i = 1:length(pList)
    
    fName = [root,'/',behRoot,num2str(pList(i)),'_PerceptualGroupingSequential.mat'];
    load(fName)
    
    allsubs.condition(i,:) = stim.condition(:)';
    allsubs.setSize(i,:) = stim.setSize(:)';
    allsubs.change(i,:) = stim.change(:)';
    
    allsubs.accuracy(i,:) = stim.accuracy(:)';    
    allsubs.rt(i,:) = stim.rt(:)';
    
    allsubs.subjects(i,:) = pList(i);
         
end

save('allsubs_data.mat','allsubs','-v7.3');







