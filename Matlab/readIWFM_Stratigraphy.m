function STRAT = readIWFM_Stratigraphy(varargin)
% STRAT = readIWFM_Stratigraphy(filename, Nnd, Nheader)
%
% Reads the stratigraphy file of IWFM
%   Nnd are the numebr of nodes
%   Nheader are the lines to skip. The value for C2VSim V1 is 105
% IMPORTANT this works only if the aquifer has 4 layers

Nnd = 30179;
Nheader = 105;
Nlay = 4;
if length(varargin) > 0
    filename = varargin{1};
end
if length(varargin) > 1
    Nnd = varargin{2};
end
if length(varargin) > 2
    Nheader = varargin{3};
end
if length(varargin) > 3
    Nlay = varargin{4};
end
frmt = ['%f %f' repmat(' %f',1,2*Nlay)];
fid = fopen(filename,'r');
strat = textscan(fid, frmt, Nnd, 'HeaderLines', Nheader);
fclose(fid);
STRAT = cell2mat(strat);
end

