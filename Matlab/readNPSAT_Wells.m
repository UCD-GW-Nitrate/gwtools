function wells = readNPSAT_Wells(filename)
% wells = readWells(filename)
% Reads the NPSAT well input file. 
% The output is a table where the rows are:
% X, Y are the well coordinates
% Top, Bot are the top and bottom elevation of the well screen
% Q is the pumping rate


fid = fopen(filename,'r');
Nwells = fscanf(fid, '%d',1);
temp = fscanf(fid, '%f',Nwells*5);
fclose(fid);
wells = reshape(temp, 5, Nwells)';
wells = [(1:size(wells,1))' wells];

wells = array2table(wells, 'VariableNames', {'Eid','X','Y','Top','Bot', 'Q'});
