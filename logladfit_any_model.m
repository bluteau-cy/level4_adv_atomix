function [epsilon,slope_ladLog,epsi_ci,epsi_ladNL,slope_ladNL,epsi_lsqNL,slope_lsqNL]=logladfit_any_model(k,Pk, modelEpsi,logslope,veldir,tNL)
% Replaces the mlefit_any_model by finding slope and/or offset of
% log-transformed data. Created to test alternate fitting methods.
% It presumes a shape: (Pk=b*k^logslope) and uses least-squa6re deviation
% rather than least-squares. The best epsilon is returned using the
% prescribed logslop, while the slope is returned by allowing both
% intercept and slope to vary in the lad.
% Required inputs:
%       k: for wavenumber (or cpm) or even frequency in Hz (Rad/s). this must be consistent with the modelSpec units 
%	P= observed spectra power in same units as modelSpec (typically in (units timeserie)^2/(rad/m)
%	logslope=-5/3 i.e. expected slope with respect "k" data
%                    Pk=a*k^logslope, which beomes log(Pk)=logslope*log(k)+a
%       modelEpsi: function/calc that returns epsilon from the calculated b (or offset in log-transformed data)
%               For example. for the vertical velocity , cte depends on the vel direction used (F_dir=[1 4/3 4/3]);    
%                       modelEpsi=@(bt)(((10.^bt)./cte).^1.5); %                           
%                       cte=1.5*(18/55)*F_dir; (F_dir=1 or 4/3 for longi and trans/vertical
%       tNL=1/0 of wherher nonlinear lad is done. 0 will speed up code
% Output:
%       epsilon: obtained from taking the median of log-transformed "lad" data
%       epsi_ladNL: obtained doing the minimzation on non-transformed data
%       slope: best slope of the spectra i.e., we now check if it's close
%          to -5/3
%__Created by CBluteau, used in getting epsilon from ADVs using inertial subrange model (Bluteau et al, 2011 L&O: methods) and also when fitting Nasmyth to shear spectrum (useful for high epsilon). See Bluteau et al 2016 JTECH for more details

%modelEpsi=@(bt)(((10.^bt)./alp).^1.5);
if nargin<6
   tNL=0; 
end
%% Log transform
ind=find(k>0);
Pk=Pk(ind);
k=k(ind);
x=log10(k); % rad/m
y=log10(Pk);


%% Lin transformed version
if veldir==1
    alp=1.5*(18/55);
else
    alp=1.5*(18/55)*4/3;
end
fitModel=@(epsi,ko)(inertial_epsi_model(ko,epsi,veldir));
ladModel=@(epsi)(1e5*sum(abs(Pk-fitModel(epsi,k))));
ladSlopeModel=@(b)(1e5*sum(abs(Pk-alp.*b(1).^(2/3).*k.^b(2))));
nlsqSlopeModel=@(b)(1e5*sum((Pk-alp.*b(1).^(2/3).*k.^b(2)).^2));
nlsqModel=@(epsi)(1e5*sum((Pk-fitModel(epsi,k)).^2)); % nl Lsq min
%% This would require generalisation
dat=y-logslope.*x;
blsq=median(dat); 
yfit=logslope.*x+blsq; % best fit using least-absolute deviation with one unknown. See my fitting paper for this regression type.
resid = y - yfit; % y-(logslope*x+blsq)= y-logslope*x-blsq. yfit+res=y
if exist('bootstrp','builtin')
tmp=bootstrp(1000, @(bootr)median(yfit+bootr-logslope.*x), resid); % bootstrapping blsq median(-logslope.*x+ newY) with newY=yfit+resid
brange=prctile(tmp,[2.5 97.5]);
epsi_ci=modelEpsi(brange);%([1 3])
else
epsi_ci=[NaN NaN];
end

% I think this is the same as doing lad for just the intercept 
%logModel=@(b)(sum(abs(y-(logslope.*x+b(1))))); 
%[blsq]=fminsearch(logModel,median(dat));

logModel=@(b)(sum(abs(y-(b(1).*x+b(2))))); % allowing both to be free
[bA]=fminsearch(logModel,[logslope blsq]);
slope_ladLog=bA(1);

epsilon=modelEpsi(blsq);
%% LAD NL
if tNL==1
    [epsi_ladNL]=fminsearch(ladModel,epsilon);
    [tmp]=fminsearch(ladSlopeModel,[epsilon logslope]);
    slope_ladNL=tmp(2);

    % LSQ NL 
    [epsi_lsqNL]=fminsearch(nlsqModel,epsilon);
    [tmp]=fminsearch(nlsqSlopeModel,[epsilon logslope]);
    slope_lsqNL=tmp(2);
else
    epsi_ladNL=NaN;
    slope_ladNL=NaN;
    epsi_lsqNL=NaN;
    slope_lsqNL=NaN;
end
