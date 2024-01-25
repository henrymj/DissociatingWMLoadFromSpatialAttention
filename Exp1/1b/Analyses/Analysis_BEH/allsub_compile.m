pList = [1:4,7:12,17:18,20:21,23,33,34,36:38];

root = pwd;
behRoot =  'Data/';

for i = 1:length(pList)
    
    fName = [root,'/',behRoot,num2str(pList(i)),'_PerceptualGroupingSequentialBalanced.mat'];
    load(fName)
    
    allsubs.condition(i,:) = stim.condition(:)';
    allsubs.setSize(i,:) = stim.setSize(:)';
    allsubs.change(i,:) = stim.change(:)';
    
    allsubs.accuracy(i,:) = stim.accuracy(:)';    
    allsubs.rt(i,:) = stim.rt(:)';
    
    allsubs.subjects(i,:) = pList(i);
         
end

save('allsubs_data.mat','allsubs','-v7.3');







