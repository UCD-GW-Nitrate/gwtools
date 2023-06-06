function BC = calcIWFM_BaryCenters(ND,MSH)
%BC = calcIWFM_BaryCenters(ND,MSH) calculates the element barycenters
%   ND is a matrix [Nnd x 3] where the columns correspond to [X Y ID]
%   MSH is a matrix [Nel x >=5] where the columds correspond to 
%   [ID ND1 ND2 ND3 ND4] any additional column is ignored

% find the positions of nodes ids in the ND variable

BC = nan(size(MSH,1),3);
for ii = 1:size(MSH,1)
    BC(ii,1) = MSH(ii,1);
    n = 4;
    if MSH(ii,5) == 0
        n = 3;
    end
    idx = MSH(ii,2:n+1);
    [~, loc] = ismember(idx, ND(:,3));
    BC(ii,2) = sum(ND(loc,1))/n;
    BC(ii,3) = sum(ND(loc,2))/n;
end