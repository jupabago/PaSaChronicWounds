inputPath = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Tiff Stacks/wtd1-02/';
outputPathGlobal = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Global Bin Tiff Stacks/wtd1-02/';
outputPathAdaptive = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Adaptive Bin Tiff Stacks/wtd1-02/';

[red,green,blue] = Stack2volume(inputPath);%split the image by the 3 channels.
[adaIm015,otsuIm, otsuHand, differ] = BinarizeAndCompare (red);
%Volume2Stack(outputPathGlobal, otsuIm);
%Volume2Stack(outputPathAdaptive, adaIm025);

function [adapI, otsuI,otsuHI, diff] = BinarizeAndCompare (volume)%This script compares the adaptive threshold with the global
adapT = adaptthresh(volume, 0.15);
adapI = imbinarize (volume, adapT);
adapI = (bwareaopen(adapI ,10));

otsuHT = CollectThresholds(volume)
otsuHI = imbinarize (volume, otsuHT);
otsuHI= (bwareaopen(otsuHI,10));

otsuI = imbinarize (volume);
otsuI = (bwareaopen(otsuI ,10));

[w,h,d] = size(volume);
matching = nnz(adapI ==otsuI);
diff = matching/(w*h*d);
end

function [redVolume, greenVolume, blueVolume] = Stack2volume(directory)
imageFolder=dir([directory '/*.tif']);%the star is for removing the two files that aren't tiffs
slices = size(imageFolder,1);
[width, height,~] = size(imread(strcat(directory,'/',imageFolder(1).name)));
[redVolume, greenVolume, blueVolume]= deal(zeros(width, height, slices));
for slice= 1:slices
    imageInt = imread(strcat(directory,'/',imageFolder(slice).name));
    image = im2double(imageInt);
    redVolume(:,:,slice) = squeeze(image(:,:,1)); 
    greenVolume(:,:,slice) = squeeze(image(:,:,2)); 
    blueVolume(:,:,slice) = squeeze(image(:,:,3));
end
end

function Volume2Stack(directory, volume)
[~,~] = mkdir(directory);
[~,~,slices] = size(volume);
for slice= 1:slices
imageName = strcat(directory,'/red_',GetSlice(slice),'.tif');%create image name to store
imwrite(volume(:,:,slice),imageName);%save image
end
end

function slice = GetSlice(idx)
if(idx>=100)
    slice =num2str(idx);
elseif(idx>=10)
    slice = strcat('0', num2str(idx));
else
    slice = strcat('00', num2str(idx));
end
end

function threshold = CollectThresholds(volume)
[~,~,slices]=size(volume);
totalCounts = 0;
for slice = 1:slices
    [counts] = imhist(volume(:,:,slice));%get histogram for this slice of the volume
    totalCounts= totalCounts+counts;     
end
threshold=otsuthresh(totalCounts);%this correction is because matlab starts at 1 and timepoints at 0
end
