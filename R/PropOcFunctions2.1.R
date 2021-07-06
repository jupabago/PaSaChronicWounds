#' Read image files from folder
#'
#' \code{sum} returns a list with the names of the image files in the folder.
#'
#'these functions were optimized for capped spheres on 10/14/19
#'
SlicesTo3D<-function(filesList,channel,side){
  image1<-readTiff(filesList[1])#read first image to get its dimensions
  image3D <- array(0, c(image1@size[1], image1@size[2], length(filesList)))#preallocate array for image given using image dimensions
  for(i in 1:length(filesList)){#keep in mind that it is not zero-index
    #print(paste("loading ch1...", i, "/", length(filesList), filesList[i]))
    image3D[,,i]<-pixmap::getChannels(readTiff(filesList[i]), colors = channel)
  }
  return(image3D)
}

IdVoxel<-function(thresholdedImage){#this makes a df with the coordinates from all voxel objects
  dfVoxelCoords <- data.frame(which(thresholdedImage == 1, T))
  colnames(dfVoxelCoords) <- c("x", "y", "z")
  return(dfVoxelCoords)
}

SampleVoxels<-function(voxelsCoords, sampleSize){#randomly samples non-empty voxels from a thresholded image
  if(dim(voxelsCoords)[1]<sampleSize){sampleSize <- dim(voxelsCoords)[1]}
  sampledVoxels <- sample(1:dim(voxelsCoords)[1], size = sampleSize, replace = T)
  sampleVoxelsCoords <- voxelsCoords[sampledVoxels,]
  return(sampleVoxelsCoords)
}

CorrectSubsetImage<-function(focus,thresholdedImage, xydim, zdim){#output is corrected outbox along with an list of 6 digits... 
  if ((focus$x-xydim)<1) lowBndX = 1 else lowBndX = focus$x-xydim #this prevents the code from going out of bounds in the lower end
  if ((focus$y-xydim)<1) lowBndY = 1 else lowBndY = focus$y-xydim
  if ((focus$z- zdim)<1) lowBndZ = 1 else lowBndZ = focus$z- zdim
  if ((focus$x+xydim)>dim(thresholdedImage)[1]) upBndX = dim(thresholdedImage)[1] else upBndX = focus$x+xydim #this prevents grabbing voxels out of bounds in high end.
  if ((focus$y+xydim)>dim(thresholdedImage)[2]) upBndY = dim(thresholdedImage)[2] else upBndY = focus$y+xydim
  if ((focus$z+ zdim)>dim(thresholdedImage)[3]) upBndZ = dim(thresholdedImage)[3] else upBndZ = focus$z+ zdim
  
  updatedAxes<-list(lowBndX=focus$x-lowBndX,upBndX=upBndX-focus$x, lowBndY=focus$y-lowBndY,upBndY=upBndY-focus$y,lowBndZ=focus$z-lowBndZ,upBndZ=upBndZ-focus$z)
  
  outbox<- thresholdedImage[(lowBndX):(upBndX), (lowBndY):(upBndY), (lowBndZ):(upBndZ)]
  #I'll leave this here for future reference...
  #print(paste("focus point: ", focus$x,focus$y,focus$z))
  #print(paste("dimensions: ",(focus$x-lowBndX),(upBndX-focus$x), (focus$y-lowBndY),(upBndY-focus$y),(focus$z-lowBndZ),(upBndZ-focus$z)))
  #print(paste("desired coordinates: ",(focus$x-xydim),(focus$x+xydim), (focus$y-xydim),(focus$y+xydim),(focus$z-zdim),(focus$z+zdim)))
  #print(paste("actual coordinates: ",(lowBndX),(upBndX), (lowBndY),(upBndY), (lowBndZ),(upBndZ)))
  return(list(Axes=updatedAxes,SubsetImage=outbox))
}

DoubleCappedSphereVolume<-function(r,h1, h2){#this also works with negative height, im not sure why #CRAZYMATHS
  fullVolume<-volume2<-4/3*pi*(r^3)
  cap1Volume <- pi*h1**2/3*(3*r-h1)
  cap2Volume <- pi*h2**2/3*(3*r-h2)
  return(fullVolume-cap2Volume-cap1Volume)#since they can't meet at this point, it doesn't matter
}

HollowDoubleCappedSphereVolume<-function(innerR,outerR,cappedR1,cappedR2){#since you can get a negative height, lets prevent that
  if(cappedR1>outerR) cappedR1<-outerR
  if(cappedR2>outerR) cappedR2<-outerR
  outerVolume<-DoubleCappedSphereVolume(outerR,outerR-cappedR1,outerR-cappedR2)
  if(cappedR1>innerR) cappedR1<-innerR
  if(cappedR2>innerR) cappedR2<-innerR

  innerVolume<-DoubleCappedSphereVolume(innerR,innerR-cappedR1,innerR-cappedR2)

  return(outerVolume-innerVolume)
}

SpherePropOcc<-function(binsList, frequencyList, voxelSize,cappedR1,cappedR2){#input are bins, 
  propOcc<-array(0,length(frequencyList))#allocate memory
  for (i in 1:length(frequencyList-1)){
    propOcc[i]<-frequencyList[i]*voxelSize/HollowDoubleCappedSphereVolume(binsList[i],binsList[i+1], cappedR1,cappedR2)#this is the prop oc
  }
  return(propOcc)
}

DistanceRadius<-function(focus,thresholdedImage, xydim, zdim, xyscale, zscale){
  #takes a voxel and collects all the surrounding voxels within the specified sphere
  #input of xydim and zdim is in voxels
  #maximum distance allowed is capped at xydim*xyscale
  #function to correct them should be here and the output is "outbox"
  correctedSubsample<-CorrectSubsetImage(focus,thresholdedImage, xydim, zdim)#this is where I correct the subsampling of the image and get the axes for the total volume
  #first convert the axes from pixels to microns cubed:
  scalingFactors<-c(xyscale,xyscale,xyscale,xyscale,zscale,zscale)#this list is to be used to correct the scale of the axes 
  #print(paste("scalingFactors:", scalingFactors))
  scaledAxes<-map2(correctedSubsample$Axes,scalingFactors, ~.x*.y)#these are the corrected search dimensions in microns instead of voxels
  outbox<-correctedSubsample$SubsetImage#this line calls the function to correctly subset the image
  outboxCoords<-(IdVoxel(outbox))#find non-empty voxels in box
  outboxCoords<-outboxCoords-1#move to the origin, since it is 1-indexed
  outboxCoords$x<-outboxCoords$x*xyscale#convert pixels to real lengths
  outboxCoords$y<-outboxCoords$y*xyscale
  outboxCoords$z<-outboxCoords$z*zscale
  #radius<-(((xydim*2)+1)/2)*xyscale#the plus one is because the total diameter of the circle is xydim+xydim+1(focus voxel)
  radius<-xydim*xyscale #this is a test run
  #under this new pipeline, this may not be the focus 
  scaledFocusPoint<-c(scaledAxes$lowBndX+xyscale, scaledAxes$lowBndY+xyscale, scaledAxes$lowBndZ+zscale)
  #print(paste("FocusPoint: ", focus))
  #print(paste("scaledFocusPoint1: ", c(focus[1]*xyscale,focus[2]*xyscale,focus[3]*zscale)))
  #print(paste("scaledFocusPoint: ", scaledFocusPoint))
  distances<-as.data.frame(t(cdist(t(as.matrix(scaledFocusPoint)), as.matrix(outboxCoords))))
  if(nrow(distances)>0){#this prevents code from breaking if there are no pixels in vecinity
    distances<-subset(distances,!(distances[1]==0))#this removes the distance to self
  }
  if(nrow(distances)>0){
    distances<-subset(distances,!(distances[1]>(radius)))#this filters out all the distances larger than the radius
  }
  #code crashes if you dont do independent "if statements"
  return(list(Distances=distances,Axes= scaledAxes))
}

RawBinningNorm<-function(dfDistances, binSize, voxelVol, listAxes){#takes the list of pairwise distances between focus and object in box and outputs frequency density table binned by distances of "binSize" in microns
  numDistances<-apply(dfDistances,1,as.numeric)#turn input df into numeric so hist() works
  bins<-seq(0, max(dfDistances)+binSize,binSize)#establish the binning criteria
  histogram<-hist(numDistances, breaks = bins, plot = FALSE)
  histCounts<-as.data.frame(histogram$counts)#these are the raw counts per bin
  #histDensity<-as.data.frame(histogram$counts/length(numDistances))#these are the proportion of counts over non-empty pixels
  histDensityFunction<-as.data.frame(histogram$density*binSize)#these are the same as the top but using the function from hist class
  #This is where the interesting part starts
  normalRadius<-listAxes[which.max(listAxes)]#this is the longest radius
  cappedRadius1<- listAxes[which.min(listAxes)]#I picked the lowest radius
  listAxes[which.min(listAxes)]<-NULL#then I remove it from the list
  cappedRadius2<-listAxes[which.min(listAxes)]#then i pick the second lowest radius
  totalVolume<-DoubleCappedSphereVolume(normalRadius[[1]],normalRadius[[1]]-cappedRadius1[[1]],normalRadius[[1]]-cappedRadius2[[1]])
  histProportion<-as.data.frame(histogram$counts*voxelVol/totalVolume)#these are the proportion of counts over all pixels in volume
  histNormSphereBinning<-as.data.frame(SpherePropOcc(bins,histogram$counts,voxelVol,cappedRadius1[[1]],cappedRadius2[[1]]))#then I add them to the function and go from there
  #adding both radii is fine because if they are the same as the search radius it wont be affected, and if they are not, then it's taken care of...
  finalData<-cbind(histogram$mids, histCounts, histDensityFunction, histProportion, histNormSphereBinning)
  colnames(finalData)<- c("Distance", "Counts", "BinDensity", "TotalDensity", "PropOc")
  return(finalData)
}

NormRawMergeStats<-function(dataList){
  allData<-bind_rows(dataList)
  allData%>%group_by(Distance)%>%summarise(MeanCounts=mean(Counts), MeanBinDensity=mean(BinDensity),MeanTotalDensity=mean(TotalDensity), MeanPropOcc=mean(PropOc), SDPropOcc = sd(PropOc) )%>% as.data.frame()->finalData
  return(finalData)
}

NormRawLoopVoxels<-function(voxelsList, thresholdedImage, xyBoxDim, zBoxDim, xyScale, zScale, loopBinSize ){
  loopBox<- vector(mode = "list", length = length(voxelsList[[1]]))#pre-allocate space
  voxelVolum<-xyScale*xyScale*zScale
  for (i in 1:length(voxelsList[[1]])){
    distancesAndAxes<-DistanceRadius(voxelsList[i,],thresholdedImage,xyBoxDim,zBoxDim,xyScale,zScale)#this just calculates the distances and gives axes based on focus
    distanceList<-distancesAndAxes$Distances#this includes just the distances
    axesList<-distancesAndAxes$Axes#this includes just the distances
    if(nrow(distanceList)>0){#this prevents code from breaking if there are no pixels in vecinity
      loopBox[[i]]<-RawBinningNorm(distanceList, loopBinSize, voxelVolum, axesList)
      #erased the sphere volume because I have to calculate it inside here anyways...
    }
  }
  return(loopBox)
}

#Pipeline
NormRawPipeline<-function(redImg,greenImg,xydimSearch,zdimSearch,sampleSize,xyRealDim,zRealDim, pipeBinSize){
  print("Identifying Voxels")
  voxCoordsRed<-IdVoxel(redImg)#this gives a list of the non-empty voxels
  voxCoordsGreen<-IdVoxel(greenImg)
  print("Sampling Voxels")
  sampleVoxRed<-SampleVoxels(voxCoordsRed, sampleSize)
  sampleVoxGreen<-SampleVoxels(voxCoordsGreen, sampleSize)
  print("Calculating distances self-self 1")
  distListRR<-NormRawLoopVoxels(sampleVoxRed, redImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Calculating distances self-other 1")
  distListRG<-NormRawLoopVoxels(sampleVoxRed, greenImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Calculating distances self-other 2")
  distListGR<-NormRawLoopVoxels(sampleVoxGreen, redImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Calculating distances self-self 2")
  distListGG<-NormRawLoopVoxels(sampleVoxGreen, greenImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Merging data")
  finalRR<-NormRawMergeStats(distListRR)
  finalRG<-NormRawMergeStats(distListRG)
  finalGR<-NormRawMergeStats(distListGR)
  finalGG<-NormRawMergeStats(distListGG)
  
  finalData<-gdata::combine(finalGG,finalRG,finalRR,finalGR)
  return(finalData)
}

NameAndMerge<-function(dataList, name){
  allData<-bind_rows(dataList)
  allData %>% mutate(source = name) %>% as.data.frame()->finalData
  return(finalData)
}

AllDataPipeline<-function(redImg,greenImg,xydimSearch,zdimSearch,sampleSize,xyRealDim,zRealDim, pipeBinSize){
  print("Identifying Voxels")
  voxCoordsRed<-IdVoxel(redImg)#this gives a list of the non-empty voxels
  voxCoordsGreen<-IdVoxel(greenImg)
  print("Sampling Voxels")
  sampleVoxRed<-SampleVoxels(voxCoordsRed, sampleSize)
  sampleVoxGreen<-SampleVoxels(voxCoordsGreen, sampleSize)
  print("Calculating distances self-self 1")
  distListRR<-NormRawLoopVoxels(sampleVoxRed, redImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Calculating distances self-other 1")
  distListRG<-NormRawLoopVoxels(sampleVoxRed, greenImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Calculating distances self-other 2")
  distListGR<-NormRawLoopVoxels(sampleVoxGreen, redImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Calculating distances self-self 2")
  distListGG<-NormRawLoopVoxels(sampleVoxGreen, greenImg, xydimSearch, zdimSearch, xyRealDim, zRealDim, pipeBinSize)
  print("Merging data")
  finalRR<-NameAndMerge(distListRR,"finalRR")
  finalRG<-NameAndMerge(distListRG,"finalRG")
  finalGR<-NameAndMerge(distListGR,"finalGR")
  finalGG<-NameAndMerge(distListGG,"finalGG")
  
  finalData<-bind_rows(finalGG,finalRG,finalRR,finalGR)
  return(finalData)
}