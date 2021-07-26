function STRAT = readIWFM_Stratigraphy(filename, Nnd, Nheader)
% STRAT = readIWFM_Stratigraphy(filename, Nnd, Nheader)
%
% Reads the stratigraphy file of IWFM
%   Nnd are the numebr of nodes
%   Nheader are the lines to skip. The value for C2VSim V1 is 105
% IMPORTANT this works only if the aquifer has 4 layers

fid = fopen(filename,'r');
strat = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', Nnd, 'HeaderLines', Nheader);
fclose(fid);
STRAT = cell2mat(strat);
end

