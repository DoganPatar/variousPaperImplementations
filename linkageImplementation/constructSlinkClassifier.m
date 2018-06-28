function [treeMatrix]=constructSlinkClassifier(features)

[pies,lambda]=slink(features);

[distSort, I] = sort(lambda);
numOfLevel = length(lambda)-1;
tempLabel = I;
tempPies = pies;

treeMatrix = zeros([numOfLevel,5]);
temp = length(lambda) + 1;
for i = 1:numOfLevel
    treeMatrix(i,1) = tempLabel(i);
    treeMatrix(i,2) = tempPies(I(i));
    treeMatrix(i,3) = distSort(i);
    treeMatrix(i,4) = I(i);
    treeMatrix(i,5) = pies(I(i));
    
    tempPies(tempPies == treeMatrix(i,1)) = temp;
    tempPies(tempPies == treeMatrix(i,2)) = temp;    
    tempLabel(tempLabel == treeMatrix(i,1)) = temp;
    tempLabel(tempLabel == treeMatrix(i,2)) = temp;    
    temp = temp+1;
end


