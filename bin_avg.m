function [ndata]=bin_avg(xdat,data,xbins)
% To bin data that is collected unevenly in depth into an evenly spaced vector
% xbins= vector of x data to bin within.. centers the data around xbins
% xdata= unevenly (or even) spaced dim (e.g depth) where the data was collected
% Improvements:
%		Consider using histcounts.m from matlab release 2015
[m n]=size(data);
if m<n % more col than rows
    data=transpose(data);
    [m n]=size(data);
end

if m~=length(xdat)
    error('xdat vector must have the same nbre of rows or cols as data');
end
%% Computations
%dd=min(depth):dbin:max(depth); % new evenly spaced depth vector
dfx=diff(xbins);
dfx(end+1)=dfx(end);

ndata=zeros([length(xbins) n]);
for ii=1:length(xbins)
    ind=find(xdat>=xbins(ii)-dfx(ii) & xdat<xbins(ii)+dfx(ii));
    tmp=nanmean(data(ind,:));
    ndata(ii,:)=tmp;
end
%dd=dd+0.5*dbin; % centering bins
%ind=find(xdat>xbins(end));
%ndata(ii+1,:)=data(end,:);