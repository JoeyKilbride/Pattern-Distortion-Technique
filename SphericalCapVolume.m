

function [V] = SphericalCapVolume(h,Rb)
%% Calculates the volume of a droplet assuming it's a spherical cap.
%% input 1: Droplet height in m (array or float)
%% input 2: Droplet base radius in m (array or float)
%% Ouput 1: Volume in Litres (array or float) 
%%
    V = (1/6)*pi.*h.*(3*Rb.^2+h.^2);
    V=V*1000; % litres
end
