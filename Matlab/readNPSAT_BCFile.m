function BC = readNPSAT_BCFile(filename)
% BC = readNPSAT_BCFile(filename) 
% Reads the main boundary file of the NPSAT flow model
fid = fopen(filename,'r');
Nbc = cell2mat(textscan(fid,'%f',1));

for ii = 1:Nbc
    C = textscan(fid,'%s %f %s',1);
    type = C{1,1}{1};
    n = C{1,2};
    value = C{1,3}{1};
    C = textscan(fid,'%f %f',n);
    XY = [C{1,1} C{1,2}];
    BC(ii,1).Type = type;
    BC(ii,1).Value = value;
    BC(ii,1).XY = XY;
end
fclose(fid);
end