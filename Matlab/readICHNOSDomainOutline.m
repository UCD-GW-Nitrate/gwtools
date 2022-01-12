function D = readICHNOSDomainOutline(filename)
% D = readICHNOSDomainOutline(filename)
%
% Reads the domain outline input file for the ICHNOS code.
% The returned structure containts the polygon verticels and the
% orientation of the polygons. 1 means that the points inside the polygon
% are inside the domain, while 0 means that the points inside the polygon
% are outside of the domain

str = fileread(filename);
lines = regexp(str, '\r\n|\r|\n', 'split')';

idx = 1;
id_cnt = 1;
while 1
    C = textscan(lines{idx,1},'%f');
    idx = idx+1;
    nv = C{1}(1);
    orient = C{1}(2);
    pp = [];
    for ii = 1:nv
        C = textscan(lines{idx,1},'%f');
        idx = idx+1;
        pp = [pp;C{1}(1) C{1}(2)];
    end
    D(id_cnt,1).ID = id_cnt;
    D(id_cnt,1).poly = pp;
    D(id_cnt,1).Orient = orient;
    id_cnt = id_cnt + 1;
    if idx > length(lines)
        break;
    end
    if isempty(lines{idx,1})
       break; 
    end
end
