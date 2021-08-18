%{
rawTifPathBase = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Tiff Stacks New/';
binTifPathBase = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Binary Images/';
aggsFilePathBase = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Aggregate lists/'; %names are different
slicedAggsFilePathBase = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Results per slice/'; %names are different

for sampleNumber = 1:2
rawTifPath = [rawTifPathBase,'mono_d1_',GetNum(sampleNumber),'/'];
binTifPath = [binTifPathBase,'mono_d1_',GetNum(sampleNumber),'/'];
aggsFilePath = [aggsFilePathBase,'mono_d1_',GetNum(sampleNumber)];
slicedAggsFilePath = [slicedAggsFilePathBase,'mono_d1_',GetNum(sampleNumber)];
xySize = .415; %.208 or .132
zSize = .52; %.365 or .361
Tiff2Data (rawTifPath,binTifPath,aggsFilePath,slicedAggsFilePath, xySize, zSize);
end
%}

for sampleNumber = 
iterateSamples('paphz_d4_',sampleNumber,.415,.52)
end
iterateSamples('papqsL_d4_',1,.208,.52)
for sampleNumber = 3:10
iterateSamples('papqsL_d4_',sampleNumber,.208,.52)
end

for sampleNumber = 1:9
iterateSamples('pawt_d4_',sampleNumber,.208,.52)
end

function iterateSamples(sampleName, sampleNum,xy,z)

rawTifPathBase = '/run/user/937561/gvfs/smb-share:server=130.207.66.142,share=raw_data/Confocal/Carolyn/2020/Chronic wounds/Tiff Stacks New/';
binTifPathBase = '/run/user/937561/gvfs/smb-share:server=130.207.66.142,share=raw_data/Confocal/Carolyn/2020/Chronic wounds/Binary Images/';
aggsFilePathBase = '/run/user/937561/gvfs/smb-share:server=130.207.66.142,share=raw_data/Confocal/Carolyn/2020/Chronic wounds/Aggregate lists/'; %names are different
slicedAggsFilePathBase = '/run/user/937561/gvfs/smb-share:server=130.207.66.142,share=raw_data/Confocal/Carolyn/2020/Chronic wounds/Results per slice/'; %names are different
rawTifPath = [rawTifPathBase,sampleName,GetNum(sampleNum),'/'];
binTifPath = [binTifPathBase,sampleName,GetNum(sampleNum),'/'];
aggsFilePath = [aggsFilePathBase,sampleName,GetNum(sampleNum)];
slicedAggsFilePath = [slicedAggsFilePathBase,sampleName,GetNum(sampleNum)];
disp(['Sample ',sampleName,GetNum(sampleNum)])
Tiff2Data (rawTifPath,binTifPath,aggsFilePath,slicedAggsFilePath, xy, z);
end

function Tiff2Data (rawTifPath,binTifPath, aggsFilePath,slicedAggsPath, xySize, zSize)
wienerSize = ceil(5/xySize); %image filter neighborhod set at 5 microns
tic
[ch1,ch2,ch3] = Stack2volume(rawTifPath);%split the image by the 3 channels.
toc
disp('filtering and stretching channel 1')
filtered1 = imbinarize(FilterImage(ch1, wienerSize));
clear ch1
toc
disp('filtering and stretching channel 2')
filtered2 = imbinarize(FilterImage(ch2, wienerSize));
clear ch2
toc
disp('filtering and stretching channel 3')
filtered3 = imbinarize(FilterImage(ch3, wienerSize));
clear ch3
toc
disp('combining and saving')

%now remove bleeding from red into the green channel and isolated pixels and create a clean volume
%%!!!Order is really important here, it is input as red, green blue, but matlab doesnt let you do it in the actual function
[cleanRed, cleanGreen, cleanBlue] = CleanExportVolume(binTifPath, filtered1,filtered2, filtered3, slicedAggsPath);
%Use clean volumes created to get aggregate sizes
toc
disp('creating aggregate structure red');
redAggList=create3dStructure(cleanRed, xySize, zSize);
resultsfilenameRed = strcat(aggsFilePath ,'_Sa.csv');
csvwrite(resultsfilenameRed,redAggList)
toc
disp('creating aggregate structure green');
greenAggList=create3dStructure(cleanGreen, xySize, zSize);
resultsfilenameGreen = strcat(aggsFilePath ,'_Pa.csv');
csvwrite(resultsfilenameGreen,greenAggList)
toc
disp('done');
end


function filteredVolume = FilterImage(volume, wienerSize)
[width, height,slices] = size(volume);
filteredVolume = zeros(width, height, slices);
for slice= 1:slices
    stretchedImg = imadjust(volume(:,:,slice));%maximizes range of intensity
    weinerImage = wiener2(stretchedImg, [wienerSize wienerSize]);%smooths image when it's low contrast and leaves it if high contrast
    filteredVolume(:,:,slice)= weinerImage;
end
clear volume %empty original volume from memory
end

function [ch1Volume, ch2Volume, ch3Volume] = Stack2volume(directory)%takes a folder with tiffs and returns a 3d-volume/matrix
imageFolder=dir([directory '/*.tif']);%the star is for removing the two files that aren't tiffs
%slices = 10; %this is in case I need to just test this step
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

function [cleanRedVol,cleanGreenVol,cleanBlueVol]=CleanExportVolume(directory,redVolume,greenVolume,blueVolume, slicedAggsDirectory)%create clean volumes and export them as tiffs
[~,~] = mkdir(directory);
[width, height,slices] = size(redVolume);
[cleanRedVol, cleanGreenVol, cleanBlueVol]= deal(zeros(width, height, slices));
slicedAggsResults = zeros(slices*3, 3);
for slice= 1:slices
redImage = bwareaopen(redVolume(:,:,slice),10);%remove isolated voxels
cleanRedVol(:,:,slice) = redImage;
greenImage = greenVolume(:,:,slice);
greenImage = greenImage - redImage; %correct for bleeding
idx = greenImage < 0;%this identifies -1's
greenImage(idx) = 0;
greenImage = bwareaopen(greenImage,10);
cleanGreenVol(:,:,slice) = greenImage; 
blueImage = bwareaopen(blueVolume(:,:,slice),10);
cleanBlueVol(:,:,slice) = blueImage;
imwrite(redImage,strcat(directory,'/R_',GetSlice(slice),'.tiff'));
imwrite(greenImage,strcat(directory,'/G_',GetSlice(slice),'.tiff'));
imwrite(blueImage,strcat(directory,'/B_',GetSlice(slice),'.tiff'));
%sliced aggregates part,
%this could be a function, but I would be passing so many arguyments that it's not even worth it. Also, I like tidy data
aggResultRow = ((slice-1)*3)+1;
%add the respective counts
slicedAggsResults(aggResultRow,1)= nnz(redImage);
slicedAggsResults(aggResultRow+1,1)= nnz(greenImage);
slicedAggsResults(aggResultRow+2,1)= nnz(blueImage);
%add species as channels since matlab doesnt like strings in matrices
slicedAggsResults(aggResultRow:aggResultRow+2,2)= 1:3;%this is nice because I can just add all those simultaneously
%add the slice number
slicedAggsResults(aggResultRow:aggResultRow+2,3)= slice;%this is nice because I can just add all those simultaneously
end
resultsfilename = strcat(slicedAggsDirectory,'.csv');
csvwrite(resultsfilename,slicedAggsResults)
end

function aggregateSizeList = create3dStructure(cleanVolume, xySize, zSize)
structure = bwconncomp(cleanVolume,18);
totalAggregates = structure.NumObjects;
aggregateSizeList = zeros(totalAggregates,1);
for agg = 1:structure.NumObjects   
    aggregateSizeList(agg,1)= GetObjectSize(structure,agg,xySize, zSize);%volume aggregate 1
end
end

function objectSize = GetObjectSize(threeDStructure,aggNumber,xySize, zSize)
objectSize = numel(threeDStructure.PixelIdxList{aggNumber});
objectSize = objectSize*xySize*xySize*zSize;%conversion of voxels to um^3
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

function imageNumber = GetNum(idx)
if(idx>=10)
    imageNumber =num2str(idx);
else
    imageNumber = strcat('0', num2str(idx));
end
end


%deprecated functions
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