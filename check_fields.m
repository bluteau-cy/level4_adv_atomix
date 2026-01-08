function nF=check_fields(Dat,bF)
% Fct to check whether all the fields in cell array bF exists or not in struc Dat.
% Output:  a cell array of the fields that exist in structure Dat
bF=unique(bF);
ind=find(isfield(Dat,bF));
if ~isempty(ind)
    nF={bF{ind}};
else
   nF=[]; 
end

end