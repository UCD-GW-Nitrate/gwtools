function polys = readICHNOSProcPolys(filename)
%polys = readICHNOSProcPolys(filename) Reads the processor polygon file

fid = fopen(filename,'r');
c = textscan(fid,'%f',1);
Npolys = c{1};
polys(Npolys,1).ID = [];
polys(Npolys,1).Coord = [];
for ii = 1:Npolys
    c = textscan(fid,'%f %f',1);
    polys(ii,1).ID = c{1};
    Nvert = c{2};
    c = textscan(fid,'%f %f',Nvert);
    polys(ii,1).Coord = [c{1} c{2}];
end

fclose(fid);

end