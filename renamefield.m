function s=renamefield(s,fldname0,newname)
%function s=renamefield(s,fldname0,newname)
%    Change the fieldnames in the structure s
%fldname0: old name
%

if strcmp(fldname0,newname)
    disp('New name is the same, skipping')
    return;
end


if isfield(s,fldname0)
    s.(newname)=s.(fldname0);
    s=rmfield(s,fldname0);
else
    warning([fldname0, ' doesn''t exist'])
end
