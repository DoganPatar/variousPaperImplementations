function [choice] = classifyOnSlink(x,treeMatrix,features)
%x          : new sample which will be classified
%treeMatrix : constructed classifier with constructSlinkClassifier
%features   : all leaves of Slink

% Function returns corresponding cluster/class of the new x sample. 

[ind,~] = size(treeMatrix);
leafSize = ind +1;
while 1
    label1 = treeMatrix(ind,1);
    label2 = treeMatrix(ind,2);
    
    br1 = features(treeMatrix(ind,4),:);
    br2 = features(treeMatrix(ind,5),:);
    disTo1 = norm(x - br1);
    disTo2 = norm(x - br2);
    
    if disTo1 <= disTo2 
        choice = label1;
    else
        choice = label2;
    end
    
    if choice <= leafSize
        break;
    else
       ind = choice - leafSize;
    end
end