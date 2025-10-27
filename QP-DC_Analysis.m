clear all
close all

global lambda;global res;global n_m;global xSize;global ySize;
path='T:\Data\AutoRAPID\AutoRAPID\20240807_LDN_Stimulation\14-19-01\';

cd(path);
cd('field_retrieval')
alpha=0.19; %0.15 for Hb, 0.19 for other biological cells
n_m=1.337;
sampleList=dir('Field*reduced*.mat');
sampleSeq=1:length(sampleList);
for sampleNum=sampleSeq
    clearvars -except path alpha n_m sampleList sampleSeq sampleNum xSize ySize res lambda
    areaList=[];
    deformationList=[];
    dryMassList=[];
    cellList=[];
    refractiveindexList=[];
    zList=[];
    aspectRatioList=[];
    intensityList=[];
    positionList=[];
    
    croppedMap={};
    croppedIntMap={};
    contourMap={};
    eventNumber=0;
    tic
    [xSize,ySize,~]=size(retPhase);

    for iter=1:size(retPhase,3)
        % Analysis
        p=squeeze(retPhase(:,:,iter));
        amplitudeMap=retAmplitude(:,:,iter);
        level = multithresh(p,1);
        
        numeFocusFlag=0;
        if level<.3
        else
            level = 0.5;
            cellMap=(abs(p)>level); %logical table. 0/1
            cellMap=bwareafilt(cellMap,[50 2000000]);
            se=strel('disk',1);
            cellMap=imdilate(cellMap,se);
            cellMap=bwareaopen(cellMap,10);
            cellMap2=bwlabel(cellMap);
            for kk=1:max(cellMap2(:))
                try
                    %
                    cellMap3=(cellMap2==kk);
                    tempContour=bwboundaries(cellMap3);
                    tempContour=tempContour{1};
                    [k av] = convhull(tempContour(:,1),tempContour(:,2));
                    tempContour=tempContour(k,:);
                    contourLength=sum(sqrt(sum((tempContour-circshift(tempContour,[1 0])).^2,2)));
                    contourLength=contourLength*res;
                    areaVal=av.*(res.^2);
                    
                    positionVal=regionprops(cellMap3,'Centroid');
                    positionVal=positionVal(1).Centroid;
                    positionVal=round(positionVal(2:-1:1));
                    
                    if numeFocusFlag==0
                        [intensityStack,zStep]=numericalFocusing(amplitudeMap,p);
                        intensityStack=(intensityStack);
                        numeFocusFlag=1;
                    end
                    [zVal,Focused_Sample]=findingSlice(intensityStack,positionVal,cellMap3);
                    Focused_Sample=(Focused_Sample);
                    zVal=zStep(zVal);
                    amplitudeVal=abs(Focused_Sample);
                    tempPhase=unwrap2(double(angle(Focused_Sample)));
                    p=tempPhase-mode(round(tempPhase(:)*1000)/1000);
                    
                    cellMap=((p)>level); %logical table. 0/1
                    cellMap=bwareafilt(cellMap,[50 2000000]);
                    cellMap=imdilate(cellMap,se);
                    cellMap=bwareaopen(cellMap,10);
                    cellMap4=bwlabel(cellMap);
                    cellMap3=cellMap4==cellMap4(positionVal(1),positionVal(2));
                    
                    hVal=regionprops(cellMap3,'MajorAxisLength');
                    vVal=regionprops(cellMap3,'MinorAxisLength');
                    aspectRatioVal=hVal.MajorAxisLength/vVal.MinorAxisLength;
                    
                    tempContour=bwboundaries(cellMap3);
                    tempContour=tempContour{1};
                    
                    [k av] = convhull(tempContour(:,1),tempContour(:,2));
                    
                    tempContour=tempContour(k,:);
                    contourLength=sum(sqrt(sum((tempContour-circshift(tempContour,[1 0])).^2,2)));
                    contourLength=contourLength*res;
                    
                    areaVal=av.*(res.^2);
                    dryMassVal=(sum(sum(cellMap3.*p)).*(res.^2))*lambda/(2*pi*alpha);
                    X=sum(sum(cellMap3.*p)).*(res.^2);
                    deltan=X*lambda*3*sqrt(pi/areaVal.^3)/(8*pi);
                    refractiveindex=deltan+n_m;
                    if ((min(tempContour(:,2))==1)||(max(tempContour(:,2))==size(p,2)))
                    else
                        areaList=[areaList;areaVal];
                        deformationList=[deformationList;1-2*sqrt(pi.*areaVal)/contourLength];
                        dryMassList=[dryMassList;dryMassVal];
                        cellList=[cellList; [iter,kk]];
                        refractiveindexList=[refractiveindexList;refractiveindex];
                        zList=[zList;zVal];
                        aspectRatioList=[aspectRatioList;aspectRatioVal];
                        positionList=[positionList;positionVal(2)];
                        
                        [tempX tempY]=find(cellMap3>0);
                        tempSizeX=[max(1,min(tempX(:))-5),min(xSize, max(tempX(:))+5)];
                        tempSizeY=[max(1,min(tempY(:))-5),min(ySize, max(tempY(:))+5)];
                        croppedMap{size(croppedMap,2)+1}=int8(p(tempSizeX(1):tempSizeX(2),tempSizeY(1):tempSizeY(2))/2/pi*128);
                        croppedIntMap{size(croppedIntMap,2)+1}=int8(amplitudeVal(tempSizeX(1):tempSizeX(2),tempSizeY(1):tempSizeY(2))/2*128);
                        contourMap{size(contourMap,2)+1}=int8(tempContour-[tempSizeX(1)-1 tempSizeY(1)-1]);
                        
                    end
                    eventNumber=eventNumber+1;
                catch e
                end
            end
            
        end
        numeFocusFlag=0;
    end
    fileName=strcat('Analysis_',sampleList(sampleNum).name);
    try
        timeVal=toc;
        save(strcat(fileName,'.mat'),'timeVal','eventNumber','res','croppedMap','croppedIntMap','contourMap','zList','areaList','positionList','deformationList','dryMassList','cellList','refractiveindexList','aspectRatioList')
    catch
    end
end

function [intensityStack,zStep]=numericalFocusing(amplitudeMap,phaseMap)
global n_m;global res;global lambda;global xSize;global ySize;
%%
zStep=-7:0.2:7;
intensityStack=(zeros(size(phaseMap,1),size(phaseMap,2),length(zStep)));
amplitudeMap=(amplitudeMap);
phaseMap=(phaseMap);
for zz=1:length(zStep)
    I_sample=(amplitudeMap.*exp(1i.*phaseMap));
    I_sample(1:21,:)=1;
    I_sample(end-21:end,:)=1;
    
    FT_I_H = fftshift(fft2((I_sample)));
    S = PSF(xSize, ySize, res, lambda,n_m, zStep(zz));
    intensityStack(:,:,zz)= ifft2(ifftshift(FT_I_H.*S));
end
end

function [zVal,Focused_Sample]=findingSlice(intensityStack,positionVal,pmap)
pmap=(imdilate(pmap,strel('disk',10)));
focusList=(zeros(size(intensityStack,3),1));
for zz=1:size(intensityStack,3)
    focusList(zz)=std(std(gradient(abs(intensityStack(:,:,zz))).*pmap));
end
[~,ind]=min(focusList);
zVal=(ind);
%%
Focused_Sample=intensityStack(:,:,zVal);
end