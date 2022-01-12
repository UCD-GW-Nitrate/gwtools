function PARAM = readIWFM_Node_Param(varargin)
% PARAM = readIWFM_Node_Param(varargin)
% or
% PARAM = readIWFM_Node_Param(filename, Nlines, Nnodes, NLayers)
% Reads the parameters of the section: 
%
% OPTION 2 - SPECIFY PARAMETERS FOR EACH MODEL NODE
% of the GROUNDWATER COMPONENT MAIN FILE
%
% The arguments are
%
% filename (Required)
% Number of line where the data start default is 54958
% Number of nodes
% Number of layers
if nargin == 1
    filename = varargin{1};
    idx = 54958;
    Nnd = 30179;
    Nlay = 4;
elseif nargin == 2
    filename = varargin{1};
    idx =  varargin{2};
    Nnd = 30179;
    Nlay = 4;
elseif nargin == 3
    filename = varargin{1};
    idx =  varargin{2};
    Nnd = varargin{3};
    Nlay = 4;
elseif nargin == 4
    filename = varargin{1};
    idx =  varargin{2};
    Nnd = varargin{3};
    Nlay = varargin{4};
end

PARAM.PKH = zeros(Nnd,4);
PARAM.PS = zeros(Nnd,4);
PARAM.PN = zeros(Nnd,4);
PARAM.PV = zeros(Nnd,4);
PARAM.PL = zeros(Nnd,4);
% read the entire file
str = fileread(filename);
% split it into lines
lines = regexp(str, '\r\n|\r|\n', 'split')';
for ii = 1:Nnd
    for j = 1:Nlay
        C = strsplit(strtrim(lines{idx,1}));
        idx = idx + 1;
        if j == 1
            PARAM.PKH(ii,j) = str2double(C{1,2});
            PARAM.PS(ii,j) = str2double(C{1,3});
            PARAM.PN(ii,j) = str2double(C{1,4});
            PARAM.PV(ii,j) = str2double(C{1,5});
            PARAM.PL(ii,j) = str2double(C{1,6});
        else
            PARAM.PKH(ii,j) = str2double(C{1,1});
            PARAM.PS(ii,j) = str2double(C{1,2});
            PARAM.PN(ii,j) = str2double(C{1,3});
            PARAM.PV(ii,j) = str2double(C{1,4});
            PARAM.PL(ii,j) = str2double(C{1,5});
        end
    end
end