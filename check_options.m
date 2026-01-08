function Opt=check_options(DefOpt,Opt)
%function Opt=check_options(mfields,def_val,Opt)
%	Created to assign default values to "Opt" when some aren't specified.
%		Was tired of repeating the code
%	Required inputs:
%		DefOpt: structure where field represent default values (entry) to apply in Opt if not specified
%			e.g., DefOpt.nu=1.36e-6;

%mfields={'intens_diff','corr','pctgd','ignorebeam'};
mfields=fieldnames(DefOpt);
%def_val=[50 100 50 0]; % defaults values corresponding to mfileds
if isempty(Opt) % use default qaqc screening parameters
    Opt=DefOpt;
else
    qafields=fieldnames(Opt);
    for jj=1:length(mfields)
        if ~any(strcmp(mfields{jj},qafields)) % mfields jj wasn't specified by the user
            warning(['You havent specified field ',mfields{jj},' so using default value of: ', num2str(DefOpt.(mfields{jj}))]);
            Opt.(mfields{jj})=DefOpt.(mfields{jj});
        end
        
        if isempty(Opt.(mfields{jj}))
            Opt.(mfields{jj})=DefOpt.(mfields{jj});
            warning([mfields{jj},' was empty so using default value of: ', num2str(DefOpt.(mfields{jj}))]);
        end
    end
end

end