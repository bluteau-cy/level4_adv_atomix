function [epsilon,  misfit, kli, allmisfit,ax,Opt]= idm_atomix_adv_spectra(v1_mean, Pw ,fw, veldir,dof,visc,Sopt)
% function [epsilon,  misfit, best_ind, allmisfit,ax,Opt]= idm_search_adv_spectra(v1_mean, Pw ,fw, veldir,dof,visc,Sopt)
% Rehacked code in 2022 to ressemble chi and shear probe code. Namely do an
% initial sweep to identify the inertial subrange.
%   Major changes:
%       - There's no extending of the fit anymore. It does one sweep on a band-avg spectra, then does
% another high-res sweep in the region (2x dec) of the identified location.
%       - MAD is the only criteria considered to find the most probable inertial subrange, provided the k median is not beyon kn~0.1.
%       - There's no rejection based on nonlinear power fit slope. This qty is calculated from the final best fit and returned though for further post-processing
%       - kLo is also no longer used to exclude wavenumbers (?)
% Required inputs
% EXPLAIN
%	fw and Pw in rad/s
%	veldir: number from 1-3 for which velocity direction to compute the spectra
%Optional
%  SOpt: structure with many fields
%	 - dec: 1 value for minimum number of decades desired for fitting. Making this small will make the code slower.
%		Accepts fractional values ie 0.5 for half a decade( default is 0.80).
%		In theory, you would prefer to fit over a decade i.e dec=1.
%	- tplot: 1/0 logical value whether or not a plot is desired. By default yes
%	- visc: numerical value for water kinematic viscosity for calculating the theoretical linit of the inertial subrange (based on epsilon and viscosity)
%			By default a value 1e-6 m2/s is used.
%   'NS': shear or stratification in rad/s for low wavenumber limit (3/Le where Le=sqrt(e/N^3) or S^3)

% Outputs
%   kli: wavenumber limits of best fit [rad/m]

% Thigs to do, mar 2022
%   Initialse ALL relevant misfit
%% Old help fix
% It must also receive the mean velocity output
% 	Spectral fitting to intertial subrange of velocity spectra using maximum likelihood method
% The best (lowest) MAD misfit portion of the spectra is retained provided
% it doesn't fill in too much into viscous subrange
%
% Portions of the spectra  that are excluded are those which yield
%	1) the median kn>0.1 (more than half of the spectra sits in the dissipation range) or
%	2) are outside the user specified frequency range to fit (optional argument)
%	3) [not implemented]median KL<3 if the stratification or shear is provided (optional)
%	5) the power fitted slope by non-linear leastsquares is within [-2.66 -.66]; Should be -5/3=-1.66
% Output:
%	 epsilon=dissipation estimate  and associated misfit critiria for the "best fit" as a structure.
%		look at epsiQAQCcalcs for the misfit criteria calculated.
%	misfit.var and misfit.mad misfit criteria as defined by Ruddick et al , 2000 (Var-Eq.23 and MAD2 Eq.24)
%	kli: wavenumbers (rad/m) of best fit.
%	ax: represents the axis handles of the best fit plot
%   Opt: options used whcih will reflect the defaults when none were
%   provided through Sopt (options)

% Improvements and Limitations
% Renamedfrom idm_SpectraSearch and simplified options in 2022 for ATOMIX,
% and mimic chi sweep routine

%% User defined stuff and defaults

DefOpt.tplot=0; % yes ploting
DefOpt.fastR=1.05; %runs faster version of  idm_SpectraSearch (will shift each fit by a few points). BEst to do more band-avg of spectra and use tfast=0
%DefOpt.visc=1e-6; % viscosity
DefOpt.dec=0.7;  % 10^dec= decadal range

DefOpt.lowfreqlim=[]; DefOpt.highfreqlim=[];  % lowest and highest freq [rad/s] to fit irrespective of theory. HArd-coded
DefOpt.knlim=0.1; % rad/m end of inertial subrange
DefOpt.kNy=100*2*pi; % max "fixed" wavenumber in rad/m
DefOpt.epsiSearch=[1e-10 1e-2];
DefOpt.kmin=0;

%% Obsolete options.
DefOpt.minR2=0; % min R2 based on non-linear least square regression
DefOpt.NS=[]; % shear/strat in rad/s
DefOpt.brange=0.75*[-1 +1]-5/3; % substract 5/3 and this yields the range of acceptable spectral slope from nonlinear power fit.


%% Check optuibs
if nargin<6
    Sopt=[];
end

if isempty(Sopt)
    Opt=DefOpt;
else
    Opt=check_options(DefOpt,Sopt);
end

%disp(log10(Opt.epsiSearch))
switch veldir
    case{1} % same as advection direction
        F_dir=1;
    case{2,3}
        F_dir=4/3;
    otherwise
        error('Invalid input, there are only 3 possible velocity directions');
end
alp=1.5*(18/55)*F_dir;
%% Main script


% Initialization
epsilon=NaN;
kli=[NaN NaN];
ax=[]; % plotting purposes only
misfit=init_struct();

%vF={'mad','slope_min','slope','slope_max'};;
if isempty(fw) || any(isnan(fw));disp('Garbage spectral estimate');return;end % Nothing in spectra


%% Data conversions to remove the 0 pt in some spectra calcs
ind=find(fw>0);
fwall=fw(ind);% in rad/s
Pwall=Pw(ind); %

% % Keep only frequencies within the low and highfreqlim
% tmpf=fw./(2*pi); % in Hz to compare with freq limit
if isempty(Opt.highfreqlim) % useful to exclude very noisy frequencies
    highfreqlim=max(fwall);
else
    highfreqlim=Opt.highfreqlim;
end
if isempty(Opt.lowfreqlim) % Assign the noisef floor to a value greater than the maximum frequency. If none was specified by the user
    lowfreqlim=fwall(1);
else
    lowfreqlim=Opt.lowfreqlim;
end
% end
ind=find(fw>=lowfreqlim & fw<=highfreqlim); % exclude noise floor and lwo freq spectra not within range
SpecToFit.P=Pw(ind);
SpecToFit.f=fw(ind); % rad/s

%%
SpecRadm.k=SpecToFit.f./v1_mean; % rad/m
SpecRadm.Pk=SpecToFit.P.*v1_mean;
SpecRadm.dof=dof;
%% First do band avg and sweep?
dT=2*pi./SpecToFit.f(1);
%epsiSearch=[1e-10 1e-3];
%% Required for log fitting

fitLogModel.modelEpsi=@(bt)(((10.^bt)./alp).^1.5); %

fitLogModel.logslope=-5/3;
fitLogModel.veldir=veldir;
fitModel=@(epsi,k)(inertial_epsi_model(k,epsi,veldir)); % for MLE

[SmoothSpec.k,SmoothSpec.Pk]=band_avg(dT,SpecRadm.k,SpecRadm.Pk) ;
SmoothSpec.dof=dof*5; % dummy var for sweeping
[epsilon,kli,allmisfit,misfit]=process_sweep(SmoothSpec,fitModel,fitLogModel,visc,Opt,1.5,Opt.epsiSearch);
%[epsiT, kmed,allmisfit]=sweep_idm_any_model(SmoothSpec,fitModel,Opt.dec,[],epsiSearch);

%% SHould I do another restricted search?
fact=10.^Opt.dec;

kRange=[kli(1)./(0.75*fact) kli(2)*0.75*fact]; % new range to sweep (less than 2x dec)
if isnan(epsilon)   
    misfit=renamefield(misfit,'MAD2','mad');
    allmisfit=renamefield(allmisfit,'MAD2','mad');
    return;
end
indR=find(SpecRadm.k>=kRange(1) & SpecRadm.k<=kRange(end));
SpecR.k=SpecRadm.k(indR,:); % rad/m
SpecR.Pk=SpecRadm.Pk(indR,:);
SpecR.dof=SpecRadm.dof;

% For now this uses epsi_lad and MAD to find the best k of the inertial subrange.
[epsilon,kli,allmisfit,misfit]=process_sweep(SpecR,fitModel,fitLogModel,visc,Opt,Opt.fastR,Opt.epsiSearch);
allmisfit=renamefield(allmisfit,'MAD2','mad');
allmisfit.fw=allmisfit.kmed.*v1_mean;
allmisfit.mad_sqrt_dof=allmisfit.mad.*sqrt(dof);
misfit=renamefield(misfit,'MAD2','mad');
misfit.mad_sqrt_dof=misfit.mad.*sqrt(dof);
%% Now plotting and geting final misfit calcs
best_ind=find(SpecR.k>=kli(1) & SpecR.k<=kli(2) & SpecR.k>0);
 misfit.npts=0;
if ~isempty(best_ind) && length(best_ind)>2
    % [epsilon, misfit]=epsiQAQCcalcs(SpecFit.f(best_ind),SpecFit.P(best_ind),v1_mean,F_dir,dof,visc,Opt.NS); % this does all the epsilon and other interesting calcs
    %misfit.npts=length(best_ind);
    misfit.fw=median(SpecToFit.f(best_ind)); % rad/s
    [epsilon,misfit.slope_lad,misfit.epsi_lad_ci]=logladfit_any_model(SpecR.k(best_ind),SpecR.Pk(best_ind), ...
    fitLogModel.modelEpsi,fitLogModel.logslope,veldir);
         
    % Plotting===================================================================
    
    if Opt.tplot
       mytitle=['\epsilon =',num2str(epsilon,'%2.1e'),', best-fit slope: ',num2str(misfit.slope_lad*3,'%2.1f'),'/3 [',...
            num2str(misfit.slope_lad_ci(1)*3,'%2.1f') ,'/3 to ',num2str(misfit.slope_lad_ci(2)*3,'%2.1f') ,'/3], velocity: ',...
            num2str(veldir)];
          %ax=idm_plot(fwall, Pwall,best_ind,v1_mean,epsilon,allmisfit,misfit,mytitle,visc,F_dir); % Plotting routine
        [ax]=idm_search_plot(fwall./v1_mean,Pwall.*v1_mean,best_ind,v1_mean,epsilon,allmisfit,mytitle,visc,veldir);
    end
    
end


end
% My subfunctions to alleviate the qaqc checks within the forward and backward search loop========


%%
function misfit=init_struct()

vF={'mad','slope_log_ci','slope_log','fw','epsi_log','epsilon'};
for jj=1:length(vF)
    misfit.(vF{jj})=NaN;
end
end

%% Band avg sweep
function [sK,sPk,nAvg]=band_avg(dT,k,Pk)
% 1./dT is freq interval (i.e. delta*f)
% k and Pk in rad/m
% Output
%   sK and sPk (smoothed spectra)
nDec=1.5; % which first dec to determine point resolution (2 decades is still pretty resolved)
tdT=64; %targeted tDT
% f=f./(2*pi);
% ind=find(f>0);
% dT=1./f(ind(1)); % 512 s?

nAvg=round(dT./tdT);
if nAvg<1
    nAvg=5;
end

dK=median(diff(k));
%kmin=k(round(nAvg/2));

ind=find(k<k(1)*10.^nDec); % want to preserve same n points per dec as first 2
decF=log10(max(k)./min(k));% decades available
if decF>nDec
    npts=1.2*decF*length(ind)./nDec;
else
    npts=length(k)./nAvg;
end



sK(:,1)=logspace(log10(dK/5),log10(k(end)),npts);
%sK(:,1)=dK:nAvg*dK:k(end);
[sPk]=bin_avg(k,Pk,sK);

%[sPk, sk]= bandavg(Pk,k,nAvg,dfb,[1 2 3],13.33,4,0); %last 3
end

% %% function powerlaw slope using Clauset et al
% function clauset(k,Pk)
% plaw=@(a,b)(a.*k.^b)
% [f, gof] = fit(k, Pk, 'power1');
% % using mle
% end
