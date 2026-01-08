function TrimName=strip_ind(NetCDF,ind)
%  TrimName=strip_ind(NetCDF,ind)
% Recovers subset of all fieldnames according to supplied ind, along the
% longest dimension
vF=fieldnames(NetCDF);
if ~isempty(ind)
    for ii=1:length(vF)
        m=size(NetCDF.(vF{ii}));
        if m(1)>=ind(end)
            switch length(m)
                case{1,2}
                     TrimName.(vF{ii})=NetCDF.(vF{ii})(ind,:);
                case{3}
                    TrimName.(vF{ii})=NetCDF.(vF{ii})(ind,:,:);
            end
        else
            TrimName.(vF{ii})=NetCDF.(vF{ii});
        end

    end  
else
    TrimName=struct([]);
end
end