function [epsilon,kli,allmisfit,misfit]=process_sweep(Spec,fitModel,fitLogModel,visc,Opt,fastR,epsiSearch) % modify me if you want to change the criteria for acception a fit
%function [epsilon,kli,allmisfit,misfit]=process_sweep(Spec,fitModel,visc,Opt,fastR,epsiSearch) % modify me if you want to change the criteria for acception a fit
% Copied subfct from idm_search_adv_spectra to enable using the tool for
%   shear probes. Plz check sweep_idm_any_model or that toolbox for inputs
% Spec is in rad/m units (or cpm)
% kli: is range of inertial subrane in rad/m (or cpm) i.e. same units as
% Spec, fitModel
% fitmodel: function handle of the form (epsi,k) e.g.,  fitModel=@(epsi,k)(inertial_epsi_model(k,epsi,veldir));
% epsiSearch: range of epsilon to search. Set to [] to use defaults.
%        Sometimes the MLE struggles with very wide ranges to search
% fastR: 1.05 how fast to sweep across
% Opt: structure with a bunch of required options
%       dec: decade to search. 0.7 to 1 are good choices
%       knlim: 0.1  max non-dimen k*eta to  fit.  
%           Can use something larger for shear probe if you want to fit the viscous 
%
%You can supply Spec, kli and knlim in cpm or rad/m. It must match though the units of fitModel 

DefOpt.kNy=100; % presumes cpm
DefOpt.kmin=0; 
DefOpt.knlim=0.1/(2*pi); % presumes cpm 
Opt=check_options(DefOpt,Opt);

%%
%fact=10.^Opt.dec;
[epsiT, kmed,allmisfit]=sweep_idm_any_model(Spec,fitModel,fitLogModel,Opt.dec,Opt.kNy,fastR,epsiSearch);
allmisfit.kmed(:,1)=kmed;
allmisfit.epsilon(:,1)=epsiT;

eta=etaK(epsiT,visc);
tt=1:1:length(epsiT);

if isfield(Opt,'knlim')
    ind=find(kmed.*eta<=Opt.knlim); %find all acceptable fits with kmed range
    tt=intersect(ind,tt);
end

ind=find(allmisfit.kli(:,1)>=Opt.kmin & allmisfit.kli(:,2)<=Opt.kNy & allmisfit.npts>=8);
tt=intersect(ind,tt);
%%
vF=fieldnames(allmisfit);
if isempty(tt) % there no acceptable fit within the prescribed k-limits
    epsilon=NaN;
    kli=[NaN NaN];
    disp('++++++ No acceptable fit +++++++')
    for ii=1:length(vF)
           misfit.(vF{ii})=NaN;
    end
    return;
end

[~,indT]=min(allmisfit.MAD2(tt));
%[~,indT]=min(abs(allmisfit.slope_lad(tt)-fitLogModel.logslope)); % this
%was biased
ind=tt(indT);

epsilon=epsiT(ind);
%kli=[kmed(ind)./(log10(2)*fact) kmed(ind)*log10(2)*fact]; % range of fitted K
kli=allmisfit.kli(ind,:);


%% COuld add here a fit calculation from Opt.kmin to kli(2)?
misfit=strip_ind(allmisfit,ind);

end