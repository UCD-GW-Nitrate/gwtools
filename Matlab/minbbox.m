function f = minbbox(x,pp,tp)
% calculate the barycenter of the points
bcx = mean(pp(:,1));
bcy = mean(pp(:,2));

% translate the points so that the origin is the barycenter
pp_t = bsxfun(@minus, pp, [bcx bcy]);
M = [cosd(x) sind(x); -sind(x) cosd(x)];
pp_r = nan(size(pp_t));
for kk = 1:size(pp_t,1)
    pp_r(kk,:) = M*pp_t(kk,:)';
end
bbx_r = [min(pp_r(:,1)) min(pp_r(:,2));...
       min(pp_r(:,1)) max(pp_r(:,2));...
       max(pp_r(:,1)) max(pp_r(:,2));...
       max(pp_r(:,1)) min(pp_r(:,2))];

Minv = inv(M);
bbx_t = nan(size(bbx_r));
for kk = 1:size(bbx_r,1)
    bbx_t(kk,:) = Minv*bbx_r(kk,:)';
end
bbx = bsxfun(@plus, bbx_t, [bcx bcy]);
if tp == 0
    f = polyarea(bbx(:,1),bbx(:,2));
else
    f = bbx;
end

