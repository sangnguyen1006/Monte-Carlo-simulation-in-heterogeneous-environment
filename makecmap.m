function cmap = makecmap(Nt)
% function cmap = makecmap(Nt)
%   Currently, set for makecmap(Nt), where Nt <= 7.
%   You can add more colors as you add tissues to makeTissueList.m.

cmap = zeros(64,3);

dj = 0.05;
for i=1:64
    j = round((i-dj)/64*(Nt-1));
    if      j<=1-dj, cmap(i,:) = [1 1 1];                    % air
    elseif  j<=2-dj, cmap(i,:) = [236/255 214/255 124/255];  % normal tissue
    elseif  j<=3-dj, cmap(i,:) = [178/255 79/255 14/255];    % bruised tissue
    end
end
