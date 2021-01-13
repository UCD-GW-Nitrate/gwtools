function STRAT = readIWFM_Stratigraphy(filename, Nnd, Nheader)
%readC2Vsim_Stratigraphy Reads the stratigraphy file of IWFM
%   Nnd are the numebr of nodes
% Nheader are the lines to skip.
% IMPORTANT this workds only if the aquifer has 4 layers

fid = fopen(filename,'r');
strat = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', Nnd, 'HeaderLines', Nheader);
fclose(fid);
STRAT = cell2mat(strat);
end

