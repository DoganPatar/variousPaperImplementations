
I = imread('cameraman.tif');
noised = imnoise(I,'gaussian',0,0.005);
rampDenoised = rampAnisodiff(noised,15,1/4,9);
dynDenoised = dynamicAnisodiff(noised,19,1/4,8,2,0.5);
combDenoised = combinedAnisodiff(noised,19,1/4,8,2,0.5);

pos = [180, 0];
fSize = 12;
font = 'DejaVuSans-BoldOblique' ;
I = insertText(I,pos,'Orig','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
noised = insertText(noised,pos,'Noised','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
rampDenoised = insertText(rampDenoised,pos,'Ramp','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
dynDenoised = insertText(dynDenoised,pos,'D-alpha','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);
combDenoised = insertText(combDenoised,pos,'Combined','BoxOpacity',0,'TextColor','white','FontSize',fSize,'Font',font);

out = [I,noised,rampDenoised,dynDenoised,combDenoised];
imshow(out,[]);