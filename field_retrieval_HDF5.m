clear all
close all

NA=.55; % Numerical aperture of the objective lens
lambda=0.52; % [um] Illumination wavelength
res=13.7/(90.9443); % [um] Pixel Size
resTarget=14/40;

%% Image load
path='T:\Data\AutoRAPID\AutoRAPID\20240807_LDN_Stimulation\14-19-01';
cd(path);
bglist=dir('*.hdf5');
errorMsg=[];
for sampleNum=1:length(bglist)
    tic
    try
        cd(path);
        tempImg=h5read(bglist(sampleNum).name,"/events/image");
        img=single(tempImg(:,:,1)');
        [ii,jj]=size(img);
        targetSizeX=round(ii*res/resTarget/2);
        targetSizeY=round(jj*res/resTarget/2);
        r=round(ii*res*NA/lambda); %Radius of the Fourier mask
        yr=round(jj*res*NA/lambda); %Radius of the Fourier mask
        %% Background image
        [xSize ySize]=size(img);
        Fimg = fftshift(fft2(img))/(xSize*ySize); %FFT
        [mi,mj]=find(Fimg==max(max(Fimg(:, round(jj*0.01):round(jj*0.49) )))); % Find a Fourier peak
        mi=round(mi-ii/2-1); mj=round(mj-jj/2-1);
        c1mask = ~(mk_ellipse(yr,r,jj,ii));
        c3mask = circshift(c1mask,[mi mj]).*~c1mask;
        
        Fimg = Fimg.*c3mask; % Apply Fourier mask
        Fimg = circshift(Fimg,-[round(mi) round(mj)]);
        
        c4mask=circshift((mk_ellipse(10,10,jj,ii)),[-30 74]);
        Fimg = Fimg.*c4mask;
        
        Fimg = Fimg(ii/2-targetSizeX+1:ii/2+targetSizeX,jj/2-targetSizeY+1:jj/2+targetSizeY);
        
        [sizeFimgx, sizeFimgy]=size(Fimg);
        Fimg = ifft2(ifftshift(Fimg))*(sizeFimgx*sizeFimgy);
        Fbg=Fimg;
        
        %% Sample image
        retPhase=zeros(round(2*targetSizeX),round(2*targetSizeY),size(tempImg,3),'single');
        retAmplitude=zeros(round(2*targetSizeX),round(2*targetSizeY),size(tempImg,3),'single');
        for iter=1:size(tempImg,3)
            cd(path);
            img=single(tempImg(:,:,iter)'); %Read hologram images
            
            Fimg = fftshift(fft2(img))/(ii*jj);
            Fimg = Fimg.*c3mask;
            Fimg = circshift(Fimg,-[round(mi) round(mj)]);
            Fimg = Fimg.*c4mask;
            Fimg = Fimg(ii/2-targetSizeX+1:ii/2+targetSizeX,jj/2-targetSizeY+1:jj/2+targetSizeY);
            [sizeFimgx, sizeFimgy]=size(Fimg);
            Fimg = (ifft2(ifftshift(Fimg)))*(sizeFimgx*sizeFimgy);
            Fimg=Fimg./squeeze(Fbg); %Phase compensation by no sample background
            retAmplitude(:,:,iter)=abs(Fimg);
            
            p=-((unwrap2(double(angle(Fimg))))); %Phase Unwrapping
            [p,coeffVal,temp]=phaseCompensation(p,1); %Phase compensation
            p=p-mode(round(p(:)*1000)/1000);
            
            pMask=abs(p)>0.5;
            pMask=bwareafilt(pMask,[5 10000000]);
            pMask=imdilate(pMask,strel('disk',5));
            [p,coeffVal,temp]=phaseCompensation(p,1,~pMask); %Phase compensation
            p=p-mode(round(p(:)*1000)/1000);
            retPhase(:,:,iter)=p;
        end
        retAmplitude=retAmplitude./mean(retAmplitude,3);
        
        tempImg=max(retPhase,[],3);
        pmap=tempImg>multithresh(tempImg);
        pmap=bwareafilt(pmap,[100 1000000]);
        ind=find(squeeze(sum(squeeze(sum((retPhase.*repmat(pmap,[1 1 size(retPhase,3)]))>0.5,1)),1))<20);
        retPhase=retPhase-mean(retPhase(:,:,ind),3); %for EFP mode;
        %%
        tempImg=max(retPhase,[],3);
        pmap=tempImg>multithresh(tempImg);
        pmap=bwareafilt(pmap,[100 1000000]);
        pmap=imerode(pmap,strel('disk',5));
        ind=find(squeeze(sum(squeeze(sum((retPhase.*repmat(pmap,[1 1 size(retPhase,3)]))>0.5,1)),1))<20);
        retPhase(:,:,ind)=[];
        retAmplitude(:,:,ind)=[];
        
        cd(path)
        mkdir('field_retrieval')
        cd('field_retrieval')
        figure(1)
        subplot(2,3,[1 3]),imagesc(squeeze(max(retPhase,[],3)),[-6 6]),axis image
        subplot(2,3,[4 6]),imagesc(squeeze(max(retPhase(end/2-10:end/2+10,:,:),[],1)),[-3 3]) % show phase image
        
        colormap('jet')
        fileName=strcat('Field_',bglist(sampleNum).name,'_','_reduced');
        save(strcat(fileName,'.mat'),'retAmplitude','retPhase','xSize','ySize','NA','lambda','res','bglist','ind');
        print('-dpng',strcat(fileName,'.png'))
    catch e
    end
    toc
end