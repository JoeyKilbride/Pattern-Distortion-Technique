function [r_f, dr_fMIN, dr_fMAX, error, p_count] = findObjs(CC, XY, R, MIN)
     
     if ~exist('MIN','var')
         % third parameter does not exist, so default it to something
          MIN = 0;
     end
    
    stats = regionprops(CC, 'centroid', 'area','MajorAxisLength','MinorAxisLength');

    centroids = cat(1,stats.Centroid);
    radii = mean([[stats.MajorAxisLength];[stats.MinorAxisLength]],1)/2;
    hold on
    rmMAX = find(radii ~=max(radii)); % removing max value which is a misdetect
    centroids=centroids(rmMAX, :);
    radii=radii(rmMAX).';
    

    % check none of the regions extend outside the ROI
    DX=(XY(1)-centroids(:,1));
    DY=(centroids(:,2)-XY(2));
    theta = atan2(DY,DX)-180;
    delX = radii.*cos(theta);
    delY = radii.*sin(theta);
    xf = centroids(:,1)+delX;
    yf = centroids(:,2)-delY;
    dot2centre = ceil(sqrt(abs(XY(1)-xf).^2+abs(XY(2)-yf).^2));
    Insideidx = find(dot2centre<R); % filter out particles outside ROI region. 
    centroids=centroids(Insideidx, :);
    %areas=areas(Insideidx);
    radii = radii(Insideidx);
    viscircles(centroids,radii, 'color', 'b');
    viscircles(XY,R, 'color', 'r');
    pause(0.01)
    drawnow;
    r_f=radii(radii>MIN);
    p_count=length(r_f);
    dr_fMIN=min(radii(radii>MIN));
    dr_fMAX=max(radii(radii>MIN));
    error = std(radii(radii>MIN))/sqrt(length(radii(radii>MIN)));
    
end