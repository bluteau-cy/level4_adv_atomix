function tt=get_kInd(k,kmin,kmax)
% k: vector of values to identify index b/w kmin and kmax
% find indices b/w kmin and kmax. If either is empty then the 1st (kmin) or last k obs (kmax) is used
% I kept changing my mind, so I only have to change it once
if isempty(kmin);kmin=k(1);end
if isempty(kmax);kmax=k(end);end

tt=find(k<=kmax & k>=kmin);

    while length(tt)<8 % want at least 8 points. Move to lower k (assumes well chosen bin size)
        if tt(1)>1
             tt=tt(1)-1:1:tt(end);
         else
         	if tt(end)==kmax
         		break;
         	else
             tt=tt(1):1:tt(end)+1;
             end
         end
    end

end