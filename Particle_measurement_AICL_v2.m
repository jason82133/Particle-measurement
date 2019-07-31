%
% NAME:
%               Particle_measurement_AICL_v2
%
% PURPOSE:
%               For batch analysis of particles in multiple images.
%
%               Individual images pass a rolling-ball & bpass filter for removing background and noises. Regions of interest (ROIs) are then
%               identified after thresholding.
%               The number/area/intensity/signal-to-background ratio
%               (SBR)/length of individual particles in a single image are calculated.
%
%               Input:
%               Single .tif images (no videos) in a folder and a subfolder
%
%               Output:
%               Analysed images and the ones with numbering in .png file in a ROI subfolder; a Count and Result .txt files.
%               Area is calculated in pixel length squared within a ROI.
%               Intensity is the sum of values of individual pixels.
%               SBR is the corrected intensity, revealing as the sum of the difference values between individual pixels and average value of
%               background, divided by the backround.
%               Length is a rough approach, calculated by thinning individual ROIs to an one-dimensional skeleton.
%               The more the value close to the diffraction limit, the more discrete and distorted the value is.
%
%
%               Require the script of Length_calculation_v2.m, bpass.m, and Matlab R2015b
%
%
%               Written by Jason C Sang, University of Cambridge
%               June 2016
%
%               Updated on 2017/04/08
%                       Optimised the calculating speed by ten-times by rewriting looping structures and the length calculation method.
%                       Approximately less than 10 sec for analysing one image.
%                       User-friendly interface. Now all parameters can be fine-tuned in the Parameter setting section.
%
%               Updated on 2017/04/08
%                       Bug fixed.
%                       Applied Gaussian fit for ROI finding before bw thresholding.
%
%               Updated on 2017/04/08
%                       Now pixel size can be altered in the setting. The previous setting was approximated to a pixel size of 225 nm.
%
%               Updated on 2018/11/22
%                       The threshold for b/w is now changing with the average background intensity. i.e. the th value
%                       becomes a relative threshold to the background of each image.
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clear all

% Append the SPT codes
addpath('C:\Pile Higher and Deeper (PhD)\Area for Data\Analysis methods\Jason coding\main_tracking source code\')


MainPath = uigetdir; % Select the folder for averaged images
FolderList = dir(MainPath); %% Require a subfolder (eg. averaged) under the initial folder
FolderList = FolderList(3:end);



% Change folder
for folderN = 1:length(FolderList)
    
    pathname = [MainPath '\' FolderList(folderN).name '\'];
    NameList = dir(pathname);
    NameList = NameList(3:end);
    dname = {NameList(:).name};
    
%% Parameter setting
    
figuredir = ['C:\Pile Higher and Deeper (PhD)\Area for Data\20190326_Linda p53 gel\' FolderList(folderN).name '\ROI\']; % Destination for the figure saving; change the magenta part

datadir = ['C:\Pile Higher and Deeper (PhD)\Area for Data\20190326_Linda p53 gel\' FolderList(folderN).name '\']; % Destination for the data saving; change the magenta part

localised_tracking = 0;                % Set 1 for tracking restricted area, specified by the centre (X,Y) and the area size; 0 for entire images.

        X = 198;                       % Coordinate of X
        Y = 405;                       % Coordinate of Y
        stepper = 40;                  % Size of ROI in pixels

mag = 10;                              % Magnification of images during calculation

ig_size = 512;                         % Image size in pixels

pixel_size = 235;                      % Pixel size for length calculation (nm)

RB_size = 10;                          % Rolling-ball filter size in pixels.

lnoise = 3;                            % Level of noise size in pixels; usually 3. May be set to 0 or false, in which case only the highpass "background subtraction" operation is performed.

lobject = 30;                          % Level of object size in pixels somewhat larger than a typical object; usually 30-50, Syn=50, PrP=30

th = 18;                               % Threshold for b/w conversion.

save_status = 1;                       % Set 1 for saving figures/data; 0 for not saving.

figure_visibility = 'on';

%%
    if localised_tracking==1 % Stepper parameters for localized tracking
       
        X = yc;
        Y = xc;
       
        x1 = xc-stepper;
        x2 = xc+stepper;
        y1 = yc-stepper;
        y2 = yc+stepper;
    end

    for imageN = 1:size(dname,2)
    %% Image input and processing
    
    display(['Processing ', num2str(imageN) ,'th image in the "',  FolderList(folderN).name, '" folder:   ', dname{imageN}])
    
        % Check existence of folder
        if save_status==1
            if exist(figuredir, 'dir')==0
                warning('off','MATLAB:MKDIR:DirectoryExists')
                mkdir(figuredir);
            end
        end

        fname = [pathname dname{imageN}];
        info = imfinfo(fname);
        images = numel(info);

        % Input image
        clear imDataOri
        
        for i = 1:images
            % Bin Up
            imTemp(:,:) = imread(fname, i);
            imTemp(imTemp == 0) = mean2(imTemp);

            if localised_tracking==1
                imDataOri(:,:) = imTemp(x1:x2,y1:y2); % if using stepper
            elseif localised_tracking==0
                imDataOri(:,:) = imTemp(:,:); % if NOT using stepper
            end
            clear imTemp
        end
        clear i

        list = [];
        boundList = [];
        iList = [];
        sbrList = [];
        roiList = [];

for i = 1:images

    se = strel('disk', RB_size);
    avg_background(i) = sum(sum(double(imDataOri)))/(ig_size * ig_size);
    
    imDataB = imtophat(imDataOri, se); % rolling-ball filtered image

    imDataG(:,:) = imDataOri - imDataB; % background per unit

    imDataBG = imresize(imDataG, mag,'nearest'); %10x background

    imDataF = imresize(imDataB, mag,'nearest'); % 10x rolling-ball filtered image

    imData = imresize(imDataOri, mag,'nearest'); % 10x original image

    
    imDataGauss = imgaussfilt(imData,3); % Gaussian filter for images

    imDataF(:,:)=bpass(imgaussfilt(imDataF,3),lnoise,lobject); % Filter for image thresholding; usually 30-50, Syn=50, PrP=30
    
    clear imDataB imDataG imData


    %% Identify and track spots

    %display('Thresholding images...')
        
    BW = im2bw(imDataF, avg_background(i)*th*10^-7); % Transform to BW for boundaries; set the threshold; usually 0.005-0.015, Syn=0.019, PrP=0.0055
    BW = imfill(BW,'holes');

    clear out intensity sbr image blank background

    % parameters for peak finding
    %th = 200; % threshold; usually 200
    %sz = 30; % size of the ROI; usually 30

    figure('visible',figure_visibility);
    imagesc(imDataOri); %Plot
    colormap(hot)
    colorbar
    hold on %Plot

    image = imDataGauss;
    blank = imDataBG;
    
    clear imT
    
    % Calculate intensity
    [imT roiN] = bwlabel(BW);
    boundaries = bwboundaries(BW, 'noholes');
    roi=[];
    
    Count = roiN;
    if save_status==1
        save ([datadir 'Count.txt'], 'Count', '-ascii', '-append');
    end
    
    display(['Found ', num2str(roiN), ' particles.'])
    
    for j = 1:roiN
        ind = find(imT==j);
        [m n] = ind2sub(size(BW), ind);
        roi{j,1} = [m n];

        intensity(j) = sum(sum(double(image(ind)))); % sum of intensity of pixels
        background(j) = sum(sum(double(blank(ind)))); % sum of intensity of pixels 
        sbr(j) = sum(sum(double(image(ind))./double(blank(ind)))); % sum of signal-to-background ratio of each pixels
        
        roiOri{j,1} = roi{j}/10+0.5;
        
        b = boundaries{j,1}/10+0.5;
        
        out(j,1) = mean(roiOri{j,1}(:,1));
        out(j,2) = mean(roiOri{j,1}(:,2));
        
        plot(b(:,2),b(:,1),'w'); %Plot
        %plot(out(j,2),out(j,1),'rx') %Plot
        
        clear m n b ind
    end
    clear j

    if(~isempty(roi))
        nlist = [1:length(roi)]';
        list = [list; out];
        boundList = [boundList; boundaries];
        roiList = [roiList; roi];
        iList = [iList; intensity'];
        sbrList = [sbrList; sbr'];

        results = [list(:,1:2) nlist];
    else
        results=[0 0 0 0];
    end
    clear roi out boundaries intensity sbr

     if save_status==1
         saveas(gcf, [figuredir dname{imageN}, '.png']);
     end

        figure %Plot
        imagesc(imDataOri) %Plot
        hold on %Plot

        for i = 1:max(results(:,3))
            frames = find(results(:,3)==i);
            plot(results(frames,2),results(frames,1),'w-') %Plot
            colormap(hot)
            colorbar
            text(results(frames(1),2)+1,results(frames(1),1)+1,['\color{white} ' num2str(i)], 'FontSize', 7) %Plot
        end
        clear i

        if save_status==1
            saveas(gcf, [figuredir 'Numbering ', dname{imageN} '.png']);
        end
end

    % Organise the lists
    for k = 1:roiN
        xx = list(:,1)==results(k,1);
        yy = list(:,2)==results(k,2);
        pos = xx&yy;
        posVal = find(pos==1);
        
        boundListN(k) = boundList(posVal(1));
        roiListN(k) = roiList(posVal(1));
        iListN(k) = iList(posVal(1));
        sbrListN(k) = sbrList(posVal(1));
    end
    clear posVal xx yy pos k nlist list boundList roiList iList sbrList roiN
        

        %% Particle analysis
        if results(1,1) ~=1
        clear a I l A Area Intensity Length SBR track_total trackNo
   
        display('Analysing...');
        
        track_total=max(results(:,3));
    
        % Define format of final results
        Result=zeros(track_total, 5);
        Result(1:end,1)=1:track_total;
    
       %% Get Length
        NList = [];
        lengthList = [];
        lengthListN = [];
  
           [m n] = size(imDataOri);
           R = [];
           N = [];

           for j = 1:max(results(:,3))
               R = [R; roiListN{j}]; % particle ROIs in the results in the descending order of numbering
               N = [N; ones(size(roiListN{j},1),1).*results(j,3)]; % corresponding particle numbers in the results
           end
           clear j posF
       
           [ll, skel] = Length_calculation_v2(R, m, n, mag);
        
           [label, number] = bwlabel(skel); % thinned particles has their own numbering; look for original numbers in the following loop
           
           for num = 1:number  % find the numbering of thinning data
               ind = find(label==num);
               Rind = sub2ind(size(label), R(:,1), R(:,2));
               pos = find(ismember(Rind, ind));
               Nnum = N(pos);
               findN(num,1) = mean(Nnum);
               clear pos ind Rind Nnum
           end
           clear num skel label number R N m n
        

           lengthList = [lengthList; ll];
           NList = [NList; findN];
           %frameList = [frameList; ones(size(ll,1),1)*i];
           clear ll findN

    for k = 1:max(NList)
        posN = find(NList==k);
        
        for j = 1:length(posN)
            lengthListN=[lengthListN; lengthList(posN(j))];
        end
        clear posN j
    end
    clear k lengthList
    

    for trackNo = 1:track_total
        frames = find(results(:,3)==trackNo);
        
        for i = 1:length(frames)
            [m n] = size(imDataOri);

            R = roiListN{frames(i)};
            I(i) = iListN(frames(i));
            
        %% Get area
            b = boundListN{frames(i)}/mag+0.5;
            
            xb=b(:,2);
            yb=b(:,1);
            a= polyarea(xb,yb);
            A(i)=a;

        %% Final results reformation
        Result(trackNo, 2) = A(i); % area in pixel squared
        Result(trackNo, 3) = iListN(frames(i))/100; % real int
        Result(trackNo, 4) = sbrListN(frames(i))/100; % signal-to-background ratio
        Result(trackNo, 5) = pixel_size*(lengthListN(frames(i))/10+0.5); % length; pixel size multiplied by unit length
        
        end
    end

        %% Export to subfolder of where the figures saved
        if save_status==1
            save ([datadir 'Result.txt'], 'Result', '-ascii', '-append');
            display(['Analysis finished. Data saved.']);
            
        elseif save_status==0
            display('Analysis finished. No data saved.');
            
        end
        close all
        end
        clear boundListN roiListN iListN sbrListN lengthListN NList frame imDataBG imDataG imDataF imDataGauss BW roi boundaries image blank out intensity sbr results
    end
end
display('Execution completed!')
toc