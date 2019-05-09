%This file is for testing and showing how to use functions.

%generating basic dataset
dataSet = ceil(10*rand([10,1000]));

tic
% Single linkage
[pies,lambda]=slink(dataSet);
tree = drawLinkage(pies,lambda);
toc

[classifier] = constructSlinkClassifier(dataSet);
x = 10*rand([1,1000]);
class = classifyOnSlink(x,classifier,dataSet);


% tic
% % Complete linkage 
% [pies,lambda]=clink(dataSet);
% tree = drawLinkage(pies,lambda);
% toc
