function [L4,Seg,L5]=level4_epsilon_adv(L3,SegOpt)
% Function that processes ADV Level2 data into spectral level 3
% It doesn't do any spectral cleaning for motion related contamination
% Required inputs:
%   L3: structure with level 3 data with the min fields
%       -  UVW_VEL_SPEC, FREQ and PSPD_REL
%       -  DOF  
%Optional input
% 	SegOpt: structure with desired segmenting options. Set to [] if
% 	defaults are desired. Not the full options are detailed under
% 	idm_search_spectra
%       - tfast/tplot i.e., options of idm_SpectraSearch
%       -avgF: cell array of ancilary (avg) fields to be written across from L2 
%          Set to [] if you'd like defaults defAvgF  used 
%           Note that some fields will always be transfered over (vF)        
% Outputs:
%   L4: structure with EPSI results for all 3 direc
%   Seg: outputs the actual options that were used..
% L5: are other outputs, during testing.

DefOpt.tplot=0; % no ploting
DefOpt.tfast=1; %runs faster version of  idm_SpectraSearch (will shift each fit by a few points). BEst to do more band-avg of spectra and use tfast=0
DefOpt.visc=1.3e-6; % constant visc. Note that if  KVISC or KVISC35 exists in L3, then the visc will be set (in priority) to KVISC, then KVISC35, and finally to the supplied option before the default here is used
DefOpt.kmin=0; % can supply an array
DefOpt.dec=0.7;
DefOpt.veldir=3; % dir of vel for EPSI_FINAL after rm bad data
defAvgF={'TEMP','PRES','SAL'};

vF={'TIME','BURST_NUMBER','TIME_BNDS','N_VEL_COMPONENT','N_BNDS'}; % fields to transfer over for NetCDF writing purposes

if isempty(SegOpt)
    Seg=DefOpt;
else
    Seg=check_options(DefOpt,SegOpt);
end

if isempty(Seg.avgF)
    Seg.avgF=defAvgF;
end
Seg.avgF=check_fields(L3,[Seg.avgF vF]);
avgF=Seg.avgF;
Seg=rmfield(Seg,'avgF');


%% Deal with kinematic visc
[nB,nV,nS]=size(L3.UVW_VEL_SPEC);

viscF=check_fields(L3,['KVISC', 'KVISC35']);
if isempty(viscF)
   visc=Seg.visc.*ones([nB 1]);
else
    visc=L3.(viscF{1}); % pref to KVISC;
    if length(visc)==1
        visc=visc.*ones([nB 1]);
    end
    Seg.visc=visc;
end

%% Create kmin array, so that it can get assigned to Seg at each segment (compatible with other tools)
if length(Seg.kmin)==1
    kmin=Seg.kmin.*ones([nB 1]);
else
    kmin=Seg.kmin; % assumes same length as nB.
end
    

%% Calc
L4.K_BNDS=NaN([nB nV 2]);
%L4.SPEC_SLOPE_CI=NaN([nB nV 2]);
L4.SPEC_SLOPE=NaN([nB nV]);
L4.N_FITTED=NaN([nB nV]);

eF={'lad'}; % epsi var to add/store
lF={'lad'};

[L4,etF]=initVar(L4,eF,'EPSI',[nB nV]);
[L4,ltF]=initVar(L4,lF,'SPEC_SLOPE',[nB nV]);

for ii=1:nB
    if rem(ii,5)==0
        disp(['===============Processing segment ',num2str(ii),'===================='])
    end
    dof=L3.DOF(ii);
    fw=L3.FREQ*2*pi; % rad/s
    
    Pw=squeeze(L3.UVW_VEL_SPEC(ii,:,:))./(2*pi);
    meanU=L3.PSPD_REL(ii);
    k=L3.FREQ./meanU;
    Seg.kmin=kmin(ii);
    for jj=1:nV
        [L4.EPSI(ii,jj),  misfit, kli,L5.Allmisfit{ii,jj}]= idm_atomix_adv_spectra(meanU, Pw(jj,:)' ,fw, jj,dof,visc(ii),Seg);                
        L4.MAD_SQRT_DOF(ii,jj)=misfit.mad.*sqrt(dof);
        L4.EPSI_CI(ii,jj,:)=misfit.epsi_lad_ci;
        if all(isfinite(kli))
            L4.K_BNDS(ii,jj,:)=kli./(2*pi); % cpm
                % Assign slopes
            L4.N_FITTED(ii,jj)=length(find(k>=L4.K_BNDS(ii,jj,1) & k<=L4.K_BNDS(ii,jj,2)));
            L4.SPEC_SLOPE(ii,jj)=misfit.slope_lad;
            for tt=1:length(lF)
                L4.(ltF{tt})(ii,jj)=misfit.(['slope_',lF{tt}]);
            end
             for tt=1:length(eF)
                L4.(etF{tt})(ii,jj)=misfit.(['epsi_',eF{tt}]);
            end
        end
    end
end

L4.eta=etaK(L4.EPSI,L3.KVISC);
L4=append_struc(L4,L3,avgF);
Seg.avgF=avgF;



end

function [L4, tF]=initVar(L4,eF,pref,nB)
% pref='EPSI','SPEC_SLOPE'
for ii=1:length(eF)
    tF{ii}=[pref,'_',upper(eF{ii})];
    L4.(tF{ii})=NaN(nB);
    
end
end
