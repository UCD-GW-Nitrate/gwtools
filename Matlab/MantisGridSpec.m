function gs = MantisGridSpec(site)
%gs = MantisGridSpec(site) returns a structure with the grid specification
%data for the selected study area
%   site is one of the following names
%   CentralValley
%   TuleRiver
%   Modesto
%
%   See also CalcGridRowCol
switch site
    case 'CentralValley'
        gs.cornerX = -223300;
        gs.cornerY = -344600;
        gs.cellSise = 50;
        gs.Nrows = 12863;
        gs.Ncols = 7046;
    case 'TuleRiver'
        gs.cornerX = 41800;
        gs.cornerY = -247550;
        gs.cellSise = 50;
        gs.Nrows = 1027;
        gs.Ncols = 1068;
    case 'Modesto'
        gs.cornerX = -111800;
        gs.cornerY = -83200;
        gs.cellSise = 50;
        gs.Nrows = 1477;
        gs.Ncols = 1428;
    otherwise
        gs.cornerX = nan;
        gs.cornerY = nan;
        gs.cellSise = nan;
        gs.Nrows = nan;
        gs.Ncols = nan;
end