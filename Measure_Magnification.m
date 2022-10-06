%% Pattern Distortion Technique - Measure_Magnification.m
 % - Measured dot magnification from .avi videos
 % - Dependencies : Mag2Height.m, findObjs.m, SphericalCapVolume.m, 
 %   Distinguishable_colors.m (here: https://uk.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors),
 %   MATLAB's Image Analysis Toolbox
 % - By: Joseph Kilbride - joey.kilbride2020@my.ntu.ac.uk
 % 

clear all;
close all;

dir='specify analysis directory'; % should contain top camera .avi
cd(dir)
filename = 'enter_filename.avi';
struct_name = "Analysis_"+filename(1:end-4)+".mat"; % creating output structure name

%% Theory Parameters
% nG, Refractive index of the gas
% nS, Refractive index of the substrate 
% nL, Refractive index of the liquid
% d1, Gap between pattern under substrate
% d2, Substrate thickness

v = VideoReader(filename);
if exist(struct_name, 'file')~=2
    %% Organise directories
    oldfolder=cd(pwd);
    cd .. 
    returndir = strcat(pwd);
    cd(dir);
  

    %% _______________Retrieve Video info___________________________
    frames = read(v); % Reads in file
    nframes = length(frames(1,1,1,:));
    ypixels = length(frames(:,1,1,1));
    xpixels = length(frames(1,:,1,1));
    fprintf("There are %i frames\n", nframes)
    fprintf("There are %i pixels vertically\n", ypixels)
    fprintf("There are %i pixels horizontally\n", xpixels)
    
    %% __________User Input (at runtime)__________________________________
    
    prompt = {'Enter camera framerate [fr/s]:','Enter starting frame number:','Enter final frame number',...
        'nG (Air)','nS (Substrate)','nL (Liquid)',...
        'd1(mm) air gap under substrate','d2(mm) substrate thickness','Radius of ROI (px)',...
        "max_Rb(mm)","min_Rb(mm)"};
    dlgtitle = 'Analysis Parameters';
    dims = [1 35];
    definput = {'1','1',num2str(nframes),'1','1.46','1.33','0','1','25'};
    rec_props = inputdlg(prompt,dlgtitle,dims,definput);
    framerate=str2double(rec_props{1});
    tzero=str2double(rec_props{2});
    depinning=str2double(rec_props{3});
    nG=str2double(rec_props{4});
    nS=str2double(rec_props{5});
    nL=str2double(rec_props{6});
    d1=str2double(rec_props{7});
    d2=str2double(rec_props{8});
    ROIR=str2double(rec_props{9});
    max_Rb=str2double(rec_props{10});
    min_Rb=str2double(rec_props{11});
    n2frames = (depinning-tzero)+1;
    Time = [0:n2frames-1]./data.topFR; % calculate time series for data.


    %% ___________Initialising values _______________
    absframe=0; % Absolute frame counter for loop

    figure
    hold on
    
    if size(frames,3)>1  
        img_ref=rgb2gray(frames(:,:,:,end));
        img=rgb2gray(frames(:,:,:,1)); 
    else
        img_ref=frames(:,:,end);
        img=frames(:,:,1); 
    end
    correct_roi = false();
    while correct_roi == false()
        [ROIXY] = SelectDropletCentres(img, ROIR);
        correct_roi=false();
        roi_answer = questdlg('Are you happy with ROI centre?', ...
            'Select', ...
            "Yes","No","No");
        if roi_answer == "No"
            clear ROIXY
           [ROIXY] = SelectDropletCentres(img, ROIR);
        elseif roi_answer == "Yes"
            correct_roi=true();
        else 
            return
        end
    end
   tic  
   fieldR=sqrt(abs(ROIXY(:,1)-xpixels).^2+abs(ROIXY(:,2)-ypixels).^2);
   dr_fMIN=zeros(1,size(ROIXY,1));
   dr_fMAX=zeros(1,size(ROIXY,1));
   dr_fSigma=zeros(1,size(ROIXY,1));
   rf_bar=zeros(1,size(ROIXY,1));
   M = zeros(n2frames, size(ROIXY,1));
   r_t = zeros(n2frames,length(n2frames));
   dr_tMIN = zeros(1,n2frames);
   pcnt=zeros(n2frames,n2frames);
   dr_tMAX =zeros(n2frames,n2frames);
   drSigma_t=zeros(n2frames,n2frames);
   fprintf('Iterating through Frames:\n')
   for iD = 1:size(ROIXY,1)     
        fprintf("\t Analysing Droplet %i", num2str(iD))
        BW_ref = imbinarize(imcomplement(img_ref)); % OTSU method thresholding
        CC_ref = bwconncomp(BW_ref);
        [r_f, dr_fMIN(iD), dr_fMAX(iD), dr_fSigma(iD), pcount_f] = findObjs(CC_ref, ROIXY(iD,:), ROIR);% omitting min, assuming bakgnd is perfect.
        %TODO: Add check background detection
        rf_bar(iD)=median(r_f); % or Mean, Median is more resilient to misdetections
        counter=0;
        
        if size(frames,3)>1
                img=imbinarize(imcomplement(rgb2gray(frames(:,:,:,1))));
            else
                img=imbinarize(imcomplement(frames(:,:,1)));
        end
        h=imshow(img);
        for frame=tzero:depinning % frames of interest
            
            cla(h)
            counter=counter+1;
            if size(frames,3)>1
                i_img = rgb2gray(frames(:,:,:,frame));
            else
                i_img =  frames(:,:,:,frame);
            end
            img=imbinarize(imcomplement(i_img)); % Binary image
            imshow(i_img)
            drawnow;
            hold on
            title("Droplet = "+num2str(iD)+", frame = "+num2str(frame))
            hold off
            CC = bwconncomp(img);
            
            [r, drMIN, drMAX, drSigma, pcount] = findObjs(CC, ROIXY(iD,:), ROIR, rf_bar(iD)); % minimum size is background average
            
            %% Write detected dot values
            pcnt(iD,counter) = pcount;
            if pcount>0.1 % check we have at least one dot
                r_t(iD,counter) = median(r); 
                dr_tMIN(iD,counter) = drMIN;
                dr_tMAX(iD,counter) = drMAX;
                drSigma_t(iD,counter) = drSigma;
            else % if not just write zerosd
                r_t(iD,counter) = 0; 
                dr_tMIN(iD,counter) = 0;
                dr_tMAX(iD,counter) = 0;
                drSigma_t(iD,counter) = 0;
            end
          
        end
            
        %% Store dectection values in structure
        dot_detect.r_t = r_t;
        dot_detect.pcnt = pcnt;
        dot_detect.dr_tMIN = dr_tMIN;
        dot_detect.dr_tMAX = dr_tMAX;
        dot_detect.drSigma_t = drSigma_t;

        %% Determine droplet shape
        M(:,iD) = r_t(iD,:)./rf_bar(iD);
        M(isnan(M(:,iD)),iD)=0; % remove NaN values from the arrays 
        h_max = Mag2Height(M(:,iD), d1,d2,nG,nS,nL,max_Rb);
        h_th.max_Rb.("iD"+num2str(iD)) = h_max;
        h_min = Mag2Height(M(:,iD), d1,d2,nG,nS,nL,min_Rb);
        h_th.min_Rb.("iD"+num2str(iD)) = h_min;
        h_th.error.("iD"+num2str(iD)) = abs(h_th.max_Rb.("iD"+num2str(iD))-h_th.min_Rb.("iD"+num2str(iD)))/2;
        h_max(h_max==0) = h_min(h_max==0); % remove zero values in upper bound with lower bound values
        h_th.mean.("iD"+num2str(iD)) = mean([h_th.max_Rb.("iD"+num2str(iD)),h_th.min_Rb.("iD"+num2str(iD))],2);
        V_th.mean.("iD"+num2str(iD)) = SphericalCapVolume(h_th.mean.("iD"+num2str(iD))*1e-3,mean([max_Rb,min_Rb])*1e-3);
        %V_th.dV.("iD"+num2str(iD)) = sqrt((pi/6)*h_th.mean.("iD"+num2str(iD))*1e-3*(6*mean([max_Rb,min_Rb])*1e-3+(h_th.mean.("iD"+num2str(iD))*1e-3)^2)*diff([max_Rb,min_Rb])/2);
    end
    %% Write Output Structure
    Rb.min=min_Rb;
    Rb.max=max_Rb;
    TheoryCnst.nG=nG;
    TheoryCnst.nS=nS;
    TheoryCnst.nL=nL;
    TheoryCnst.d1=d1;
    TheoryCnst.d2=d2;
    TheoryCnst.Rb=Rb;
    
    ROI.R=ROIR;
    ROI.Centre=ROIXY;
    Analysis.Time=Time;
    Analysis.M = M;
    Analysis.h_th = h_th;
    Analysis.Depinning_frame = depinning;
    Analysis.Starting_frame = tzero;
    Analysis.PmtsTheory = TheoryCnst;
    Analysis.ROI = ROI;
    Analysis.V_th=V_th;
    Analysis.dot_detect=dot_detect;
else 
    %% Analysis already exists
    % allows replotting from generated output structure
    fprintf("\nAnalysis already exists, loading previous ...\n")
    load(struct_name)
    img = read(v,1);
    if size(img,3)>1  
        img=imbinarize(imcomplement(rgb2gray(img)));
    else
        img=imbinarize(imcomplement(img));
    end
    
end
%% Plotting Data
save(struct_name,'Analysis');
clr=distinguishable_colors(size(Analysis.ROI.Centre(:,1),1));% plotting colours
imshow(img);
set(gca,'position',[0 0 1 1],'units','normalized');
hold on
for iD = 1:size(Analysis.M,2)
    viscircles(Analysis.ROI.Centre(iD,:),50,'color',clr(iD,:), 'LineWidth',8);
    text(Analysis.ROI.Centre(iD,1),Analysis.ROI.Centre(iD,2),"iD"+num2str(iD), 'Color', 'r','FontSize',20, 'FontWeight', 'bold')
end
hold off
tempinset = getframe(gcf);
insetimage = frame2im(tempinset);
imwrite(insetimage,'Locator.png')
close all
figure
hold on


edit_plot = questdlg('Do you want to plot all the data points?', ...
            'Select', ...
            "Yes","No","No");
if edit_plot == "No"
    prompt = {'Entre number of points to average over:','Entre plotting increment:'};
    dlgtitle = 'Plotting Parameters';
    dims = [1 35];
    definput = {'20','30'};
    rec_props = inputdlg(prompt,dlgtitle,dims,definput);
    n_movmean=str2double(rec_props{1});
    nth_point=str2double(rec_props{2});
else
    n_movmean=1; % no averaging
    nth_point=1; % all points
end

X_series = Analysis.Time(1:nth_point:end); % get subset of datapoints
for iD = 1:size(Analysis.M,2)  
    
    Y_series = movmean(Analysis.h_th.mean.("iD"+num2str(iD))(:), n_movmean); % average the data
    Y_series = Y_series(1:nth_point:end); % plot subset of datapoints
    delY_series = Analysis.h_th.error.("iD"+num2str(iD))(1:nth_point:end); % get subset of datapoints
    alpha(h,0.2) % error bars transparent 
    scatter(X_series(1:length(Y_series)), Y_series,...
       36,clr(iD,:),'filled','DisplayName',"Droplet: "+num2str(iD));    
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData',...
       [h.Line.ColorData(1:3); 255*0.1]);% transparent error bars
end

main_ax = gca;

xlabel("Time (s)")
ylabel("Height (mm)")
axes('pos',[0.65 0.7 0.3 0.2])
imshow(insetimage)
hold off

figure
hold on
for iD = 1:size(Analysis.M,2)  

    Y_series = movmean(Analysis.V_th.mean.("iD"+num2str(iD))(:), n_movmean)/Analysis.V_th.mean.("iD"+num2str(iD))(1); % average the data
    Y_series = Y_series(1:nth_point:end); % plot subset of datapoints    
    plot(X_series(1:length(Y_series)), Y_series,...
       'color',clr(iD,:),'DisplayName',"Droplet: "+num2str(iD));    
    
end
ylabel("V/V_o (\mu L)")
xlabel("Time (s)")
hold off




