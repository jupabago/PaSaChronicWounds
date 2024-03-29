
%the beauty of this is that I dont have to stitch
ResultsPath = '/Users/jupabago/Documents/Whiteley/PROJECTS/SaPa Chronic Wounds/Images data';
MainPath= '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/';
%cowtd1m1p1a = Pipeline(ResultsPath,MainPath,'co','wt','1', '(11-10-20)', '1','7','1','cowtd1m1p1a');
%cowtd1m1p2a = Pipeline(ResultsPath,MainPath,'co','wt','1', '(11-10-20)', '1','7','2','cowtd1m1p2a');
%cowtd1m1p3a = Pipeline(ResultsPath,MainPath,'co','wt','1', '(11-10-20)', '1','7','3','cowtd1m1p3a');

%cowtd1m2p1a = Pipeline(ResultsPath,MainPath,'co','wt','1', '(11-10-20)', '2','15','1','cowtd1m2p1a');
%cowtd1m2p2a = Pipeline(ResultsPath,MainPath,'co','wt','1', '(11-10-20)', '2','15','2','cowtd1m2p2a');
%cowtd1m2p3a = Pipeline(ResultsPath,MainPath,'co','wt','1', '(11-10-20)', '2','15','3','cowtd1m2p3a');

cophzd1m1p1a = Pipeline(ResultsPath,MainPath,'co','phz','1', '(11-9-20)', '1','11','1','cophzd1m1p1a');
%cophzd1m1p2a = Pipeline(ResultsPath,MainPath,'co','phz','1', '(11-9-20)', '1','11','2','cophzd1m1p2a');
%cophzd1m1p3a = Pipeline(ResultsPath,MainPath,'co','phz','1', '(11-9-20)', '1','11','3','cophzd1m1p3a');

%cophzd1m2p1a = Pipeline(ResultsPath,MainPath,'co','phz','1', '(11-9-20)', '2','14','1','cophzd1m2p1a');
%cophzd1m2p2a = Pipeline(ResultsPath,MainPath,'co','phz','1', '(11-9-20)', '2','14','2','cophzd1m2p2a');
%cophzd1m2p3a = Pipeline(ResultsPath,MainPath,'co','phz','1', '(11-9-20)', '2','14','3','cophzd1m2p3a');

function imagePath = Pipeline(resultsPath,mainPath, condition,strain,day,date, mouse, ImgNum,position, sampleName)
imagePath = [mainPath, condition,' ',strain,'_day',day,' ',date, '/', 'Mouse ', mouse, '/Image ',ImgNum,' Block ', position, '_Stitch'];

[eslices, bytes,Folder] = GetStackData(imagePath);
%to iterate throught the image, do this:
iterateData(eslices,imagePath, Folder,resultsPath,sampleName); 
end

function results = iterateData (slices, path, imageFolder,resultsfilePath, sampleName)
for channel =1:3 
results = zeros(slices,12);%this is here for the simple version of the code without overlapping aggregate data
for slice= 1:slices
    sliceData = GetImageData(imread(strcat(path,'/',imageFolder(slice).name)), channel);
    results(slice,:) = sliceData;
end
resultsfilename = strcat(resultsfilePath,'/',sampleName,'_',num2str(channel),'.csv');
csvwrite(resultsfilename,results)
end
end
    
function [slices, bits, imageFolder] = GetStackData(path)
imageFolder=dir([path '/*.tif']);
slices=size(imageFolder,1);
imageName = strcat(path,'/',imageFolder(1).name);
info = imfinfo(strcat(path,'/',imageFolder(1).name));%gets image Data
bits = info(1).BitDepth/info(1).SamplesPerPixel;%gets bit depth of image
end

%{
hist0 = imhist(squeeze(I0(:,:,1)));
hist2 = imhist(squeeze(I2(:,:,1)));
hist10 = imhist(squeeze(I10(:,:,1)));
hist11 = imhist(squeeze(I11(:,:,1)));

imageData2 = GetImageData(I2);
imageData0 = GetImageData(I0);
imageData10 = GetImageData(I10);
imageData11 = GetImageData(I11);
%}

%testGetImageHistogram = GetImageHistogram(Ib, 1, 1,1,3,2);
function histogram = GetImageHistogram(image, channel, slice, timePoint, position, tile )
[imCounts,imBins ] = imhist(image);
histogram = [imCounts,imBins, ones(256,1)*channel, ones(256,1)*slice, ones(256,1)*timePoint, ones(256,1)*position, ones(256,1)*tile]; 
end

function imageData = GetImageData(im, channel)
image= squeeze(im(:,:,channel));
[width, height] = size(image);
numpixels = width *height;
[imCounts,~ ] = imhist(image);%get the histogram
zeroInt = imCounts(1);%count the instances at zero intensity
zeroPercent = zeroInt/numpixels; %count the percent zeros on an image
maxInt = imCounts(size(imCounts));
[~,firstZero] = min(imCounts);
tailEnd = image>=firstZero;%this thresholds the image after the first instance of zero
numtailEnd = nnz(tailEnd);%This is the amount of pixels beyond the first zero
totalZeroes = nnz(~imCounts);
tailBins = 257-firstZero-totalZeroes; %this is the number of intensity bins after the first zero freq that have non-zero freq. 
level = graythresh(image);
calculatedG = stretchlim(image);
meanIntensity = mean2(image);
sdIntensity = std2(image);
optimalG = getOptimalG(imCounts,numpixels);

imageData = [zeroInt,zeroPercent, maxInt,firstZero,totalZeroes,level,calculatedG(2),tailBins, numtailEnd,optimalG,meanIntensity,sdIntensity]; 
end

function optimalG = getOptimalG( counts,totalpixels)
seq = 0:1/255:1;
bins = seq';
bins(1)=[];
counts(1)=[];
optimalG = -sum(counts.*(log2(bins)))/totalpixels;
end