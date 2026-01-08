function [varm, MAD2,madshear,medLog]=ruddick_misfit(Pw,Pt)
% Ruddick et al misfit criteria
% Required inputs:
%  Pw: Spectral Observations
%  Pt: The theoretical spectraum

PwPt=Pw./Pt;
meanPwPt=nanmean(PwPt);
MAD2=nanmean(abs(PwPt-meanPwPt));% Ruddick et al
varm=var(PwPt);



%logdev=(Pw-Pt)./sqrt(Pw);
logdev=log(Pw)-log(Pt);
%madshear=nanmean(abs(logdev));
%madshear=nanmean(logdev.^1);
madshear=nanmean(abs(logdev.^1));
%madshear=var(Pw./Pt);

medLog=nanmedian(abs(logdev));