function [ax]=idm_search_plot(k,Pk,best_ind,U,epsilon,misft,mytitle,visc,veldir)
% Plots various results from sweeping a spectra, as per ATOMIX tests
%   kcyc: rad/m
% 	Pxx : in rad/m (observed spectra)
% 	mytitle: is desired title for the figure, which I now place as text
%Optional
% 	misfit: 1 value (anything) if  the bottom misfit panel is not desired. Otherwise pass the misfit data in this argument.
%Output
%	h=returns the figure handles
%			h(2) refer to those on the misfit plot
%			h(1) refers to the best fit PSD plot
% Improments and Limitations
%	Script is messy but has improved during 2015.
% CBluteau fixed some issues with kn labels in 2021
%_________________________________________________________

F_dir=[1 4/3 4/3];

alp=F_dir(veldir)*1.5*18/55;

[m,n]=size(Pk);
if all([m n]>1)
  kp=repmat(k,[1 n]);
  leglab={'u','v','w'};
else
    kp=k;
    leglab{1}=['obs ',num2str(veldir)];
end
%% Main script


ax(1)=subplot(5,1,1:2);
hF=loglog(kp,Pk.*kp.^(5/3),'.-');
hold on;

yli=get(gca,'ylim');
xli=get(gca,'xlim');
xtic=get(gca,'xtick');ytic=get(gca,'ytick');
ylabel('\Phi\timesk^{5/3} [(m/s)^2 (rad/m)^{2/3}]');

if ~isnan(epsilon)
    
    k=k(best_ind); %in rad/m
    Pt= alp*(epsilon.^(2/3)).*k.^(-5/3);% rad/m
    %pp=find(misft.kmed>=median(ft),1,'first');
    [medkt,pp]=find_nearest(misft.kmed,median(k));
    %medft=misft.kmed(pp);
    
    
    hF(end+1)=loglog(k,Pt.*k.^(5/3),'color',[0 0 0],'linewidth',3);
    leglab{end+1}=['Best fit ',num2str(veldir)];
    hl=legend(hF,leglab,'location','northwest');
    
    eta=etaK(epsilon,visc); %(visc.^3./epsilon).^0.25;
    kn=eta*xtic;
    %Ekn=ytic*U./(epsilon*visc.^5)^0.25;
    %fkolm= 0.1*U./((visc.^3./epsilon).^0.25); %1./(10*sqrt(1e-6/epsilon)); % freq higher than this are in inertial subrange
    knlim=[kn(1) kn(end)];%*eta./U;
    knlim=10.^[round(log10(knlim(1))) round(log10(knlim(2)))];
    
    ax(1)=gca;
    ax(2)=axes('Position',get(ax(1),'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right','xscale','log',...
        'Color','none', 'xlim',knlim,'ytick',[]);
    set(ax(1),'xlim',xli,'xtick',xtic,'ylim',yli,'ytick',ytic,'ygrid','on','yminorgrid','off');
    %set(ax(2),'xtick',kn,'ytick',Ekn,'xticklabel',kn,'yticklabel',Ekn);
    set(get(ax(2),'xlabel'),'string','k\eta')
    %set(ax(2),'xtick',[]);
else
    pp=[];
    medkt=[];
end


text(0.05,0.1,mytitle,'units','normalized')
%title(mytitle);
%set(ax,'xlim',xli,'ylim',yli);

ax(3)=subplot(5,1,3);
semilogy(misft.kmed,misft.epsilon ,'.-');
hold on
plot(medkt,misft.epsilon(pp) ,'o')
ylabel(' \epsilon');
%legend({'epsi\_mle','epsi\_lad','epsi\_log'},'location','southwest')


ax(4)=subplot(5,1,4);
ht=plot(misft.kmed,misft.slope_lad,'-');hold on
%plot(misft.kmed,misft.slope_log_ci(:,1),'-','color',[1 0 .5]);
%ht(2)=plot(misft.kmed,misft.slope_log_ci(:,2),'-','color',[1 0 .5]);
plot(medkt,misft.slope_lad(pp),'o')
plot(xli,[0 0],'k-');plot(xli,-5/3*[1 1],'k--');
legend(ht,{'slope','slope\_ci'},'location','northwest')
ylabel('Slope')


ax(5)=subplot(5,1,5);
loglog(misft.kmed,misft.mad_sqrt_dof,'.-');
hold on
plot(medkt,misft.mad_sqrt_dof(pp),'o')
plot(xli,2.82*[1 1],'k-');
ylabel('MAD \times DOF^{1/2}')


%     ax(6)=subplot(7,1,6);
%     semilogy(misft.kmed,misft.kL,'-');
%     hold on;
%     plot(medft,misft.kL(pp),'o')
%     plot(xli,3*[1 1],'k-');
%     ylabel('kL')
%
%     ax(7)=subplot(7,1,7);
%     plot(misft.kmed,misft.R2,'-');
%     hold on;
%     plot(medft,misft.R2(pp),'o')
%     ylabel('R^2')
xlabel('k [rad/m]');

set(ax(3:end-1),'xticklabel',[]);
set(ax(3:end),'xscale','log','xtick',xtic,'xlim',xli)
set(ax(4),'ylim',[-8/3 2/3]) % -5/3 [
set(ax(3),'ygrid','on','yminorgrid','off')
set(ax(5),'ylim',[0 5]);