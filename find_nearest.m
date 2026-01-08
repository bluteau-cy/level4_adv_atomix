function [closestVal,ind]=find_nearest(data,val)
%Syntax  [vv,ind]=find_nearest(data,val)
% Finds the nearest index (ind)  (and value vv) in the "data" vector to the
% requested number "val"
[c, ind] = min(abs(data-val));
closestVal = data(ind); % Finds first one only!