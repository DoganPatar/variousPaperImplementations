function [treeMatrix]= drawLinkage(pies,lambda)

[distSort, I] = sort(lambda);
numOfLevel = length(lambda)-1;
tempLabel = I;
treeMatrix = zeros([numOfLevel,3]);
temp = length(lambda) + 1;
for i = 1:numOfLevel
    treeMatrix(i,1) = tempLabel(i);
    treeMatrix(i,2) = pies(I(i));
    treeMatrix(i,3) = distSort(i);
    
    pies(pies == treeMatrix(i,1)) = temp;
    pies(pies == treeMatrix(i,2)) = temp;    
    tempLabel(tempLabel == treeMatrix(i,1)) = temp;
    tempLabel(tempLabel == treeMatrix(i,2)) = temp;    
    temp = temp+1;
end
dendrogram(treeMatrix);

end