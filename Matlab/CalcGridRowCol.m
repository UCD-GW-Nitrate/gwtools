function [R, C] = CalcGridRowCol(P, gridSpec)
%[R, C] = CalcGridRowCol(P, gridSpec) Calculates the rows and columns that
%correspond to coordinates
%
%   P [nx2] : are the coordinates to calculate the corresponding row and
%          column
%   gridSpec is a structure with the following fields
%       cornerX   : is the X coordinate of the left lower point of the grid.
%       cornerY   : is the Y coordinate of the left lower point of the grid.
%       cellSise   : is the cell dimension
%       Nrows   : is the number of rows
%       Ncols   : is the number of columns
%
%   See also MantisGridSpec

XG = gridSpec.cornerX + [0 cumsum(gridSpec.cellSise*ones(1,gridSpec.Ncols))];
YG = gridSpec.cornerY + [0 cumsum(gridSpec.cellSise*ones(1,gridSpec.Nrows))];

Np = size(P,1);
ind = [1:Np]';

R = nan(Np,1);
C = nan(Np,1);
P_bck = P;

id = find(P(:,1) < XG(1));
if ~isempty(id)
    C(ind(id),1) = 1;
    P(id,:) = [];
    ind(id,:) = [];
end

id = find(P(:,1) > XG(end));
if ~isempty(id)
    C(ind(id),1) = 1;
    P(id,:) = [];
    ind(id,:) = [];
end

for ii = 1:length(XG) - 1
    id = find(P(:,1) >= XG(ii) &  P(:,1) <= XG(ii+1));
    if ~isempty(id)
        C(ind(id),1) = ii;
        P(id,:) = [];
        ind(id,:) = [];
    end
end

P = P_bck;
ind = [1:Np]';
clear P_bck;
Ny = length(YG) - 1;

id = find(P(:,2) < YG(1));
if ~isempty(id)
    R(ind(id),1) = Ny;
    P(id,:) = [];
    ind(id,:) = [];
end

id = find(P(:,2) > YG(end));
if ~isempty(id)
    R(ind(id),1) = 1;
    P(id,:) = [];
    ind(id,:) = [];
end

for ii = 1:Ny
    id = find(P(:,2) >= YG(ii) &  P(:,2) <= YG(ii+1));
    if ~isempty(id)
        R(ind(id),1) = Ny - ii+1;
        P(id,:) = [];
        ind(id,:) = [];
    end
end

end