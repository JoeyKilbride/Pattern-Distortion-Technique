function [ROIXY] = SelectDropletCentres(image, ROI_radius)
%% Select ROIs on image with mouse. 
%% Input: image or frame
%% Input: Radius of ROI
%% Output: Array of ROU centres 
%% Output: Array of ROI radii

    BW = imbinarize(imcomplement(image)); % OTSU method thresholding
    imshow(BW)
    another=true();
    iD=0;
    while another == true()
       iD=iD+1;% count droplet centres selected
       roi = drawpoint(); 
       ROIXY(iD,:) = roi.Position; % xmin, ymin, width, height
       ROIR(iD) = ROI_radius;
       more_droplets = questdlg('Add another?', ...
           'Select', ...
           "Yes","No","No");

       if more_droplets == "No"
           another = false();
       elseif more_droplets == 'Yes'
           disp("Please select another centre.")
       end
    end
end