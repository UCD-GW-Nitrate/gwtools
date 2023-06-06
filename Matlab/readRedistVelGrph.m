function [Vel, Graph] = readRedistVelGrph(filename)

str = fileread(filename);
% split it into lines
lines = regexp(str, '\r\n|\r|\n', 'split')';
Nl = length(lines);
Vel.XYZ = nan(Nl,3);
Vel.proc = nan(Nl,1);
Vel.diam = nan(Nl,1);
Vel.ratio = nan(Nl,1);
Vel.VXYZ = nan(Nl,3);

idx = 1;
for ii = 1:Nl
    if isempty(lines{ii})
        continue;
    end
    C = strsplit(strtrim(lines{ii,1}));
    CC = cellfun(@str2double, C);
    Vel.XYZ(idx,:) = CC(1,1:3);
    Vel.proc(idx,1) = CC(1,4);
    Vel.diam(idx,1) = CC(1,5);
    Vel.ratio(idx,1) = CC(1,6);
    Vel.VXYZ(idx,:) = CC(1,7:9);
    idx = idx + 1;
end
if idx <= Nl
    Vel.XYZ(idx:end,:) = [];
    Vel.proc(idx:end,:) = [];
    Vel.diam(idx:end,:) = [];
    Vel.ratio(idx:end,:) = [];
    Vel.VXYZ(idx:end,:) = [];
end

fn = strsplit(filename, '.');
filegraph = [fn{1,1} '.grph'];

str = fileread(filegraph);
% split it into lines
lines = regexp(str, '\r\n|\r|\n', 'split')';

Nl = length(lines);
Graph.XYZ = nan(Nl,3);
Graph.NeighCells{Nl,1} = [];
Graph.VellCell{Nl,1} = [];

idx = 1;
for ii = 1:Nl
    if isempty(lines{ii})
        continue;
    end
    C = strsplit(strtrim(lines{ii,1}));
    CC = cellfun(@str2double, C);
    Graph.XYZ(idx,:) = CC(1,1:3);
    ncels = CC(1,4);
    Graph.NeighCells{idx,1} = CC(6:6+ncels-1);
    Graph.VellCell{idx,1} = CC(6+ncels:end);
    idx = idx + 1;
end
if idx <= Nl
    Graph.XYZ(idx:end,:) = [];
    Graph.NeighCells(idx:end,:) = [];
    Graph.VellCell(idx:end,:) = [];
end

