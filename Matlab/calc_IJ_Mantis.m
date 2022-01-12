function [II, JJ] = calc_IJ_Mantis( varargin )
%calc_IJ_Mantis returns the I J of the pnt. 
% The first input is a matrix Nx2 with the point coordinate in 3310
% The second input which is optional is a structure with information about
% the grid.
% Leave it empty for the default grid

if nargin == 1
    pnt = varargin{1};
    opt.X = -223300:50:129000;
    opt.Y = -344600:50:298550;
elseif nargin > 1
    pnt = varargin{1};
    opt = varargin{2};
end

Np = size(pnt,1);

ind = [1:Np]';

II = nan(Np,1);
JJ = nan(Np,1);
pnt_bck = pnt;

id = find(pnt(:,1) < opt.X(1));
if ~isempty(id)
    JJ(ind(id),1) = 1;
    pnt(id,:) = [];
    ind(id,:) = [];
end
id = find(pnt(:,1) > opt.X(end));
if ~isempty(id)
    JJ(ind(id),1) = 1;
    pnt(id,:) = [];
    ind(id,:) = [];
end
for ii = 1:length(opt.X) - 1
    id = find(pnt(:,1) >= opt.X(ii) &  pnt(:,1) <= opt.X(ii+1));
    if ~isempty(id)
        JJ(ind(id),1) = ii;
        pnt(id,:) = [];
        ind(id,:) = [];
    end
end

pnt = pnt_bck;
ind = [1:Np]';
clear pnt_bck;
Ny = length(opt.Y) - 1;

id = find(pnt(:,2) < opt.Y(1));
if ~isempty(id)
    II(ind(id),1) = Ny;
    pnt(id,:) = [];
    ind(id,:) = [];
end
id = find(pnt(:,2) > opt.Y(end));
if ~isempty(id)
    II(ind(id),1) = 1;
    pnt(id,:) = [];
    ind(id,:) = [];
end

for ii = 1:Ny
    id = find(pnt(:,2) >= opt.Y(ii) &  pnt(:,2) <= opt.Y(ii+1));
    if ~isempty(id)
        II(ind(id),1) = Ny - ii+1;
        pnt(id,:) = [];
        ind(id,:) = [];
    end
end




