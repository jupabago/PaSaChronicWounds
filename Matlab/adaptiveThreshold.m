rawTifPath = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Tiff Stacks/wtd1-02/';
binTifPath = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Binary Images/wtd1-02/';
aggsFilePath = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Aggregate lists/wtd1-02';
%
tic
[ch1,ch2,ch3] = Stack2volume(rawTifPath);%split the image by the 3 channels.
toc
disp('filtering and stretching channel 1')
filtered1 = imbinarize(FilterImage(ch1));
toc
disp('filtering and stretching channel 2')
filtered2 = imbinarize(FilterImage(ch2));
toc
disp('filtering and stretching channel 3')
filtered3 = imbinarize(FilterImage(ch3));
toc
disp('combining and saving')
sdsds

%now remove bleeding from red into the green channel and isolated pixels
%and create a clean volume
[cleanRed, cleanGreen, cleanBlue] = CleanExportVolume(binTifPath, filtered3,filtered1, filtered2);%Order is really important here
%Use clean volumes created to get aggregate sizes
toc
disp('creating aggregate structure red');

redAggList=create3dStructure(cleanRed);
toc
disp('creating aggregate structure green');
greenAggList=create3dStructure(cleanGreen);
%}
resultsfilenameRed = strcat(aggsFilePath ,'_red.csv');
csvwrite(resultsfilenameRed,redAggList)
resultsfilenameGreen = strcat(aggsFilePath ,'_green.csv');
csvwrite(resultsfilenameGreen,greenAggList)
toc
disp('done');

function filteredVolume = FilterImage(volume)
[width, height,slices] = size(volume);
filteredVolume = zeros(width, height, slices);
for slice= 1:slices
    stretchedImg = imadjust(volume(:,:,slice));
    weinerImage = wiener2(stretchedImg, [10 10]);
    filteredVolume(:,:,slice)= weinerImage;
end
end

function [ch1Volume, ch2Volume, ch3Volume] = Stack2volume(directory)%takes a folder with tiffs and returns a 3d-volume/matrix
imageFolder=dir([directory '/*.tif']);%the star is for removing the two files that aren't tiffs
%slices = 10;
slices = size(imageFolder,1);
[width, height,~] = size(imread(strcat(directory,'/',imageFolder(1).name)));%third dimension here is 3, one per channel
[ch1Volume, ch2Volume, ch3Volume]= deal(zeros(width, height, slices)); %initialize matrix with zeros
for slice= 1:slices
    imageInt = imread(strcat(directory,'/',imageFolder(slice).name));%read image. Intensity is not 0-1 yet
    image = im2double(imageInt);%convert to 0-1 intensity values so that filters work and thresholding work as expected
    ch1Volume(:,:,slice) = squeeze(image(:,:,1)); %separate channels. remember channels don't necessarily match colors in order
    ch2Volume(:,:,slice) = squeeze(image(:,:,2)); 
    ch3Volume(:,:,slice) = squeeze(image(:,:,3));
end
end

function [cleanRedVol,cleanGreenVol,cleanBlueVol]=CleanExportVolume(directory,redVolume,greenVolume,blueVolume)%create clean volumes and export them as tiffs
[~,~] = mkdir(directory);
[width, height,slices] = size(redVolume);
[cleanRedVol, cleanGreenVol, cleanBlueVol]= deal(zeros(width, height, slices));
for slice= 1:slices
redImage = bwareaopen(redVolume(:,:,slice),10);
cleanRedVol(:,:,slice) = redImage;
greenImage = greenVolume(:,:,slice);
greenImage = greenImage - redImage; %correct for bleeding
idx = greenImage < 0;%this identifies zeros
greenImage(idx) = 0;
greenImage = bwareaopen(greenImage,10);
cleanGreenVol(:,:,slice) = greenImage; 
blueImage = bwareaopen(blueVolume(:,:,slice),10);
cleanBlueVol(:,:,slice) = blueImage;
imwrite(redImage,strcat(directory,'/R_',GetSlice(slice),'.tiff'));
imwrite(greenImage,strcat(directory,'/G_',GetSlice(slice),'.tiff'));
imwrite(blueImage,strcat(directory,'/B_',GetSlice(slice),'.tiff'));
end
end

function aggregateSizeList = create3dStructure(cleanVolume)
structure = bwconncomp(cleanVolume,18);
totalAggregates = structure.NumObjects;
aggregateSizeList = zeros(totalAggregates,1);
for agg = 1:structure.NumObjects   
    aggregateSizeList(agg,1)= GetObjectSize(structure,agg);%volume aggregate 1
end
end

function objectSize = GetObjectSize(threeDStructure,aggNumber)
objectSize = numel(threeDStructure.PixelIdxList{aggNumber});
objectSize = objectSize*.415*.415*.52;
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