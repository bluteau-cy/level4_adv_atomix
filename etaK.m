function eta=etaK(epsi,visc)
% Function to calculate the Kolmogorov length scale
%eta=etaK(epsi,visc)
eta=(visc.^3./epsi).^0.25;

