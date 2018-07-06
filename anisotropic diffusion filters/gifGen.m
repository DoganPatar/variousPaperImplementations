clear all;
close all;
clc;

I = imread('cameraman.tif');
noised = imnoise(I,'gaussian',0,0.005);

pos = [180, 0];
fSize = 12;
font = 'DejaVuSans-BoldOblique' ;

Iout = insertText(I,pos,'Orig','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
noisedOut = insertText(noised,pos,'Noised','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);

numIm = 25;
step = 1;
for i = 1:numIm
    itter = i +step*(i-1);
    rampDenoised = rampAnisodiff(noised,itter,1/4,9);
    dynDenoised = dynamicAnisodiff(noised,itter,1/4,8,2,0.5);
    combDenoised = combinedAnisodiff(noised,itter,1/4,8,2,0.5);
    
    rampDenoised = insertText(rampDenoised,pos,'Ramp','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
    dynDenoised = insertText(dynDenoised,pos,'D-alpha','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
    combDenoised = insertText(combDenoised,pos,'Combined','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);    
    
    out{i} = [Iout,noisedOut,rampDenoised,dynDenoised,combDenoised];
end

filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:numIm
    [A,map] = rgb2ind(out{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.20);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.20);
    end
end