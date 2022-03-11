function STRAT = readIWFM_Stratigraphy(filename, Nnd, Nlay, Nheader)
% STRAT = readIWFM_Stratigraphy(filename, Nnd, Nlay, Nheader)
%
% Reads the stratigraphy file of IWFM
%   Nnd is the numebr of nodes
%   Nlay is the number of layers
%   Nheader are the lines to skip. The value for C2VSim V1 is 105

frmt = ['%f %f' repmat(' %f',1,2*Nlay)];
fid = fopen(filename,'r');
strat = textscan(fid, frmt, Nnd, 'HeaderLines', Nheader);
fclose(fid);
STRAT = cell2mat(strat);
end

