function [Pt,alp]=inertial_epsi_model(k,epsi,fdir)
%function Pt=inertial_epsi_model(k,epsi,fdir)
% Inertial subrange model for velocity.
% Req inputs:
%   k in rad/m
%   epsi
%   fdir: direction 1, 2, 3 for longi, transverse of vertical

if fdir>1 % trans or vertical
    aFac=4/3;
else
    aFac=1;
end

alp=1.5*(18/55)*aFac;% kte for theoretical spectrum
Pt= alp.*k.^(-5/3)*(epsi.^(2/3));