function ar = multipoly_area(X, Y)
%multipoly_area calculates the area of multipolygons
%
%   ar = multipoly_area(X, Y)
% where X and Y can be multiple polygons separated by nan similar to
% shapefiles
%
% Clockwise polygons are added and counterclockwise polygons are subtracted

[Xs, Ys] = polysplit(X, Y);
ar = 0;
for k = 1:length(Xs)
    a = polyarea(Xs{k,1}, Ys{k,1});
    if ispolycw(Xs{k,1}, Ys{k,1})
        ar = ar + a;
    else
        ar = ar - a;
    end
end
end