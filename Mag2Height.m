function [heights] = Mag2Height(MAG, d1,d2,n1,n2,n3,Rb)
% Calculates the matrix transfer theory height from a magnification value
% by solving a polynomial of the form ah^2+bh+c=0.
% Inputs: array/single magnification values/s
%           d1, object distance (mm)-  under glass slide
%           d2, substrate thickness (mm)
%           n1, refractive index (air) 
%           n2, refractive index (substrate) 
%           n3, refractive index (droplet liquid) 
%           Rb, droplet base radius (mm)
% Output: outputs corresponding heights in same form as input

    heights=[];% array to hold solutions
    for M_i = 1:length(MAG) % loop needed here otherwise roots satisfy all eqns.
        M=MAG(M_i);
        a=M*((2*n1/n3)-1)-1;
        b=2*M*(1-(n3/n1))*(d1+(d2*n1/n2));
        c=(M-1)*Rb^2;
        r=roots([a b c]);
        if r(1)<0 % if negative then other root is height
            h_value=r(2);
        elseif r(2)<0 % if negative then other root is height
            h_value = r(1);
        elseif r(1)>0 && r(2)>0
            h_value=min(r); % smallest value i
        elseif ~isreal(r(1)) && ~isreal(r(2))
            fprintf("\nError: Both roots are complex, returning 0 at %i... \n",M_i)
            h_value=0;
        else
            fprintf("\nError: Both roots are negative, returning 0 at %i... \n",M_i)
            h_value=0;
        end
        if isnan(h_value)% remove any NaN values
            fprintf("\nWarning: NaN values removed with zeros.")
            h_value=0;
        end
        if any(imag(r)>0) % remove any complex values
            h_value=0;
        end
        heights=cat(1,heights,h_value);
    end
    
   
   
end
