function [STREAMS, NODES] = readIWFM_StreamSpec(filename)
%[STREAMS, NODES] = readIWFM_StreamSpec(filename) 
% 
% Reads the specification file. This reads the 4.2 version and it may need
% modifications in both section where the stream nodes are defined 
% and in the rating table section
%
% STREAMS is a list of individual stream reaches 
% The IBDR variable seems unused
%
% NODES is the list of river nodes and their associated rating tables.




% read the entire file
str = fileread(filename);
% split it into lines
lines = regexp(str, '\r\n|\r|\n', 'split')';

% Get the number of streams and number of points in the rating table
current_line = 0;
for ii = 1:length(lines)
    if ~(strcmp(lines{ii}(1), '*') || strcmp(lines{ii}(1), 'C') || strcmp(lines{ii}(1), '#'))
        current_line = ii;
        break;
    end
end
C = strsplit(strtrim(lines{current_line,1}));
NRH = str2double(C{1,1});
current_line = current_line + 1;
C = strsplit(strtrim(lines{current_line,1}));
NRTB = str2double(C{1,1});

STREAMS(NRH,1).ID = [];
STREAMS(NRH,1).Name = [];
STREAMS(NRH,1).ND = [];
STREAMS(NRH,1).IBDR = [];

% Read the stream nodes specifications

for ii = 1:NRH
    while true
        current_line = current_line + 1;
        if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
            break;
        end
    end
    C = strsplit(strtrim(lines{current_line,1}));
    STREAMS(ii,1).ID = str2double(C{1,1});
    Nnodes = str2double(C{1,2});
    STREAMS(ii,1).IBDR = str2double(C{1,3});
    STREAMS(ii,1).Name = joinName({C{4:end}});
    ND = table('size',[Nnodes, 2],'VariableTypes',{'double','double'},'VariableNames',{'IRV','IGW'});
    % Read the nodes
    while true
        current_line = current_line + 1;
        if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
            break;
        end
    end
    for j = 1:Nnodes
        C = strsplit(strtrim(lines{current_line+j-1,1}));
        ND.IRV(j) = str2double(C{1,1});
        ND.IGW(j) = str2double(C{1,2});
    end
    STREAMS(ii,1).ND = ND;
    current_line = current_line + Nnodes - 1;
end
NriverNodes = STREAMS(NRH,1).ND.IRV(end);

% Read the Rating table conversion parameters
while true
    current_line = current_line + 1;
    if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
        break;
    end
end
% I havent found any reason to read this so I just skip
current_line = current_line + 3;

NODES(NriverNodes,1).IRV = [];
NODES(NriverNodes,1).BOTR = [];
NODES(NriverNodes,1).RTBL = [];

% Read the rating tables for each node
for ii = 1:NriverNodes
    while true
        current_line = current_line + 1;
        if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
            break;
        end
    end
    C = strsplit(strtrim(lines{current_line,1}));
    NODES(ii,1).IRV = str2double(C{1,1});
    NODES(ii,1).BOTR = str2double(C{1,2});
    if length(C) == 4
        RTBL = table('size',[NRTB, 2],'VariableTypes',{'double','double'},'VariableNames',{'HRTB','QRTB'});
        has_Wper = false;
    elseif length(C) == 5
        RTBL = table('size',[NRTB, 3],'VariableTypes',{'double','double','double'},'VariableNames',{'HRTB','QRTB','WRTB'});
        has_Wper = true;
    end
    RTBL.HRTB(1) = str2double(C{1,3});
    RTBL.QRTB(1) = str2double(C{1,4});
    if has_Wper
        RTBL.WRTB(1) = str2double(C{1,5});
    end
    for j = 2:NRTB
        current_line = current_line + 1;
        C = strsplit(strtrim(lines{current_line,1}));
        RTBL.HRTB(j) = str2double(C{1,1});
        RTBL.QRTB(j) = str2double(C{1,2});
        if has_Wper
            RTBL.WRTB(j) = str2double(C{1,3});
        end
    end
    NODES(ii,1).RTBL = RTBL;
end

end

function name = joinName(C)
    name = C{1};
    for ii = 2:length(C)
        name = [name ' ' C{ii}];
    end
end

