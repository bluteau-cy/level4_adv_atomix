function StruR=append_struc(StruR,StrutoAppend,vF,tforce)
% function StruR=append_struc(StruR,StrutoAppend,vF)
% function grabs all the fields in StrutoAppend, and appends them
% to structure StruR. It will not rewrite those that exist
% Optional
%   vF: cell array of fields to append. If noneprovided or empty, all the
%   fields in StruToAppend will be added

if nargin<3
    tforce=0; % 1 forces a rewrite of the variable 
    vF=[];
else
    if nargin<4
        tforce=0;
    end
end


if isempty(vF)
    vF=fieldnames(StrutoAppend);
end

%%
for ii=1:length(vF)
    if isfield(StruR,vF{ii}) && tforce==0
        disp(['Field ',vF{ii},' exists, so skipping'])
    else
        if isfield(StrutoAppend,vF{ii})
            StruR.(vF{ii})=StrutoAppend.(vF{ii});
        else
            disp(['Field ',vF{ii},' doesnt exists'])
        end
    end
end
