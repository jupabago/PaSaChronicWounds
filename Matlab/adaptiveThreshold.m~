mainPath= '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/';
fullPath = '/Volumes/raw_data/Confocal/Carolyn/2020/Chronic wounds/Tiff Stacks/wtd1-01/';

[red,green,blue] = Stack2volume(fullPath);
T = adaptthresh(red);


function [redVolume, greenVolume, blueVolume] = Stack2volume(directory)
%function  Stack2volume(directory)
imageFolder=dir([directory '/*.tif']);%the star is for removing the two files that aren't tiffs
slices = size(imageFolder,1);
%slices =5;
[width, height,~] = size(imread(strcat(directory,'/',imageFolder(1).name)));
[redVolume, greenVolume, blueVolume]= deal(zeros(width, height, slices));
for slice= 1:slices
    image = imread(strcat(directory,'/',imageFolder(slice).name));
    redVolume(:,:,slice) = squeeze(image(:,:,1)); 
    greenVolume(:,:,slice) = squeeze(image(:,:,2));
    blueVolume(:,:,slice) = squeeze(image(:,:,3));
end
end

function writeImageStack