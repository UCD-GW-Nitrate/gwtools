function wells = MantisReadWells(filename)
fid = fopen(filename);
tline = fgetl(fid);
CC = textscan(tline, '%f',1);
Nwells = CC{1};
ww = zeros(Nwells,8);
for ii = 1:Nwells
    tline = fgetl(fid);
    CC = textscan(tline, '%f %f %f %f %f %f %f %f');
    ww(ii,:) = cell2mat(CC);
end
fclose(fid);
wells = array2table(ww,'VariableNames',{'ID','X','Y','D','SL','Q', 'Ratio','Angle'});
end

