%generating basic dataset
dataSet = ceil(10*rand([10,1000]));

tic
%single linkage
[pies,lambda]=slink(dataSet);
tree = drawLinkage(pies,lambda);
toc

[classifier] = constructSlinkClassifier(dataSet);
x = 10*rand([1,1000]);
class = classifyOnSlink(x,classifier,dataSet);


% tic
% % %complete linkage 
% [pies,lambda]=clink(dataSet);
% tree = drawLinkage(pies,lambda);
% toc