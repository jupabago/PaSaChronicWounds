path ='/Volumes/Seagate Backup Plus Drive/Good images/';

%This is how to do it with the thresholds calculated from analyzing the
%thresholds from previous analysis
batchData(path, slices, timepoints, positions, SaThreshold5, 'resultsByTime7')

Pa14wtThreshold = ones(18,1)*0.081486;
Pa14mutThreshold = ones(18,1)*0.077379;
function batchData(path,slices, timepoints, positions, SaThresh, outfolder) 
GetImageDataByTime(path, '3-13-19', '62x_Salac_Pa14wt_SaPa14wt=1-1,1-10,100-1,10-1_co_SCFM2_tile2x2_3-13-19', slices, timepoints, positions, SaThresh,outfolder);
GetImageDataByTime(path, '3-19-19', '62x_Salac_Pa14wt_SaPa14wt=1-1,1-10,100-1,10-1_co_SCFM2_tile2x2_3-19-19', slices, timepoints, positions, SaThresh,outfolder);                                           
GetImageDataByTime(path, '4-17-19', '62x_Salac_Pa14wt_SaPa14wt1-11-10100-110-1_co_SCFM2_tile2x2_4-17-19', slices, timepoints, positions, SaThresh,outfolder);
                                           
GetImageDataByTime(path, '4-24-19', '62x_Salac_Pa14pqsLclean_SaPa14pqsLclean=1-1,1-10,100-1,10-1_co_SCFM2_tile2x2_4-24-19', slices, timepoints, positions, SaThresh,outfolder);
GetImageDataByTime(path, '4-25-19', '62x_Salac_Pa14pqsLclean_SaPa14pqsLclean=1-1,1-10,100-1,10-1_co_SCFM2_tile2x2_4-25-19', slices, timepoints, positions, SaThresh,outfolder);
GetImageDataByTime(path, '5-8-19', '62x_Salac_Pa14pqsLclean_SaPa14pqsLclean=1-1,1-10,100-1,10-1_co_SCFM2_tile2x2_5-8-19', slices, timepoints, positions, SaThresh,outfolder);
end
function GetImageDataByTime(path, date, name, slices, timepoints,positions, thR, outputFolder)
tic 
for position = 0:positions
        %[thR,thG] = CollectThresholds(path, date, name, slices, timepoints, position);
        ThresholdImage(path, date, name, slices, timepoints, position,thR,thG);
        toc
end
ThresholdImage(path, date, name, slices, timepoints, 2,thR,outputFolder);
toc
end

function GetImageData(path, date, name, slices, timepoints, positions)
tic    
    for position = 0:positions
        [thR,thG] = CollectThresholds(path, date, name, slices, timepoints, position);
        ThresholdImage(path, date, name, slices, timepoints, position,thR,thG);
        toc
    end
end



function ThresholdImage(path, date, name, slices, timepoints, position, redThreshold,outfolder)
imagesfilepath = strcat(path, date,'/imagesByTime');%declare name of directory to put images in
[~,~] = mkdir(imagesfilepath);%create directory if it hasn't been created. double wiggly thing is to prevent it from throwing warning
resultsfilePath = strcat(path, date,'/',outfolder);
[~,~] = mkdir(resultsfilePath);

for timepoint = 0:timepoints
    
    for slice = 0:slices
        filename = strcat(path, date,'/', name, '/', name, '_z', GetSlice(slice), '_t', GetSlice(timepoint),'_p', num2str(position));            
        I = stitchImage(filename);%collects tile of 4 images and opens in a range of 0 to 1
        %extract individual channels...
        %single channel
        ImR = squeeze(I(:,:,1));
        %ImG = squeeze(I(:,:,2));
        %ImB = squeeze(I(:,:,3));
        %binarize using "global" algorithm and global threshold
        ImRiB = imbinarize(ImR,redThreshold(timepoint+1));
        %ImGiB = imbinarize(ImG,greenThreshold(timepoint+1));
        %clean up image using maximum object size
        ImNeR = (bwareaopen(ImRiB,10));
        %ImNeG = (bwareaopen(ImGiB,10));
        singleColorResults(slice+1,1)= nnz(ImNeR);%this puts the red pixels in the results spreadsheet
        %singleColorResults(slice+1,2)= nnz(ImNeG);%this puts the green pixels in the results spreadsheet
        singleColorResults(slice+1,3)= slice;%this puts the slice number in the results spreadsheet
        %{ 
        %this part is to save the images, but I don't want to do this at this point
        %combine channels
        rgbImG = cat(3,ImNeR,ImNeG,ImB);
        imageName = strcat(imagesfilepath,'/t',GetSlice(timepoint),'_p',num2str(position),'_s',GetSlice(slice),'.tif');%create image name to store
        imwrite(rgbImG,imageName);%save image
        %}   
    end
    resultsfilename = strcat(resultsfilePath,'/t',GetSlice(timepoint),'_p',num2str(position),'.csv');
    csvwrite(resultsfilename,singleColorResults)
end
end
