function [epsiT, kmed,misfit]=sweep_idm_any_model(SpecRadm,fitmyModel,fitLogModel,dec,kNy,fastR,chiSearch)
%function [epsiT, kmed,misfit]=sweep_idm_any_model(SpecRadm,fitModel,dec,kNy,chiSearch)
%Meant to zip/search through the entire spectra and track the chit (or epsilon) as we change k of the fit.
% Required inputs:
%   SpecRadm is structure with fields:
%		 -Pk and k in rad/m
%		 -dof: degrees of freedom
%	dec: decade size over which to do the fit
%   	kNy: max k to search *really for speed, set it to max(k) if you don't want to do that
%	fitModel= function handle for model to fit in the form
%		@(chiT,k)(modelTurbSpec(epsi,chiT,k));or  @(epsi,k)(modelTurbSpec(fDir,epsi,k)) with fDir: constant?;
%       where it is called as	Pt=fitModel(chiT,k);
%
% Optional
%   kNy= max wavenumber to search set to [] or max K in the SpecRadm to
%   skip this
%   fastR: factor to move the window during search, If fastR=1.05,  log10(1.05)=0.02 so move search  window by  1/50th of decade)
%   chiSearch: is  initial "Search" domain for the MLE (default is [1e-12
%   1e-2])
% Improvements and Limitations
%   Perhaps initialize outputs
%if nargin<7
%fastR=1.1;%; %1.05; % log10(1.05)=0.02 so move search  window by  1/50th of decade)
%end
% TO-DO
% All extra variables, have the words "epsi" instead of chi. These should
% be changed
if isempty(kNy)
    kNy=max(SpecRadm.k);
end

if nargin<6
    chiSearch=[1e-12 1e-2]; % search zone
end

if isempty(chiSearch)
    chiSearch=[1e-12 1e-2]; % search zone
end
%%

fact=10^dec; % 3.166 decade to search
tfast=1;
%maxk=150; % 150 cpm is nyquist

%%
ind=find(SpecRadm.k>0);
k=SpecRadm.k(ind);
Pk=SpecRadm.Pk(ind);

Nb=length(k);

cc=1;
mn=find(k>fact*k(cc),1);
if isempty(mn);mn=Nb;end

tt=0;
while ~isempty(mn)  
    tt=tt+1;
    
    ind=get_kInd(k,k(cc),k(mn));
    misfit.npts(tt,1)=length(ind);
    
    kmed(tt)=median(k(ind));
    misfit.npts=length(ind);
    misfit.kli(tt,1:2)=[k(ind(1)) k(ind(end))];
    SpecObs.k=k(ind);
    SpecObs.P=Pk(ind);
    tmpModel=@(chiS)(fitmyModel(chiS,SpecObs.k));
    %[chiT(tt),misfit.err(tt)]=mle_any_model(SpecObs,SpecRadm.dof,chiSearch, tmpModel, 0);
%     [misfit.epsi_mle(tt,1),misfit.err(tt,1)]=mle_any_model(SpecObs,SpecRadm.dof,chiSearch, tmpModel, 0);
%     
%     [misfit.epsi_log(tt,1), misfit.epsi_log_ci(tt,:),...
%         misfit.slope_log(tt,:),misfit.slope_log_ci(tt,:),]=logfit_any_model(SpecObs.k,SpecObs.P, fitLogModel.modelEpsi,fitLogModel.logslope);
    
    [misfit.epsi_lad(tt,1),misfit.slope_lad(tt,:),misfit.epsi_lad_ci(tt,:)]=logladfit_any_model(SpecObs.k,SpecObs.P, ...
    fitLogModel.modelEpsi,fitLogModel.logslope,fitLogModel.veldir);
 
   % misfit.epsi_medianlog(tt,:)=misfit.epsi_log_ci(tt,2);
    
    epsiT(tt)=misfit.epsi_lad(tt,1);
    if isnan(epsiT(tt))
        misfit.var(tt,1)=NaN;
        misfit.MAD2(tt,1)=NaN;
    else
        Pt=fitmyModel(epsiT(tt),SpecObs.k);
        [misfit.var(tt,1), misfit.MAD2(tt,1),misfit.MAD_rolf(tt,1),misfit.MAD_other(tt,1)]=ruddick_misfit(Pk(ind),Pt);
    end
    
    if kmed(tt)>kNy
        break;
    end
    
    cc=cc+1;
    if tfast
        if k(cc)<=fastR*k(cc-1)
            cc=find(k>fastR*k(cc-1),1,'first');
        end
    end
    mn=find(k>fact*k(cc),1);
end


disp([num2str(cc),' swept portions'])
ind=find(epsiT<chiSearch(1));
epsiT(ind)=NaN;
ind=find(epsiT>chiSearch(2));
epsiT(ind)=NaN;

% ind=find(misfit.npts<8); % less than 8 pts in the fitting
% epsiT(ind)=NaN;