function [DIV_PROP,DIV_GRP] = readIWFM_DivSpec(filename)
%[DIV_PROP,DIV_GRP] = readIWFM_DivSpec(filename)
%   Reads the diversions specification file


% read the entire file
str = fileread(filename);
% split it into lines
lines = regexp(str, '\r\n|\r|\n', 'split')';

current_line = 0;
for ii = 1:length(lines)
    if ~(strcmp(lines{ii}(1), '*') || strcmp(lines{ii}(1), 'C') || strcmp(lines{ii}(1), '#'))
        current_line = ii;
        break;
    end
end
C = strsplit(strtrim(lines{current_line,1}));
NRDV = str2double(C{1,1});

Vtypes = {'double','double','double','double','double','double','double',...
          'double','double','double','double','double','double','double',...
          'string','string'};
Vnames = {'ID','IRDV','ICDVMAX','FDVMAX','ICOLRL','FRACRL','ICOLNL',...
          'FRACNL','TYPDSTDL','DSTDL','ICOLDL','FRACDL','ICFSIRIG','ICADJ',...
          'NAME','Notes'};

DIV_PROP = table('Size',[NRDV 16],'VariableTypes',Vtypes,'VariableNames',Vnames);

for ii = 1:NRDV
    while true
        current_line = current_line + 1;
        if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
            break;
        end
    end

    C = strsplit(strtrim(lines{current_line,1}));
    numeric_values = cellfun(@str2double,C);
    DIV_PROP(ii,1:14) = mat2cell(numeric_values(1:14),1,ones(1,14));
    if length(C) > 14
        DIV_PROP.NAME(ii) = C{15};
    end
    if length(C) > 15
        if length(C)>16
            DIV_PROP.Notes(ii) = C{17};
        else
            DIV_PROP.Notes(ii) = C{16};
        end
    end
end

while true
    current_line = current_line + 1;
    if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
        break;
    end
end
C = strsplit(strtrim(lines{current_line,1}));
NGRP = str2double(C{1,1});
DIV_GRP(NGRP,1).ID = [];
DIV_GRP(NGRP,1).NELEM = [];
DIV_GRP(NGRP,1).IELEM = [];
for ii = 1:NGRP
    while true
        current_line = current_line + 1;
        if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
            break;
        end
    end
    C = strsplit(strtrim(lines{current_line,1}));
    ID = str2double(C{1,1});
    NELEM = str2double(C{1,2});
    IELEM = nan(NELEM,1);
    IELEM(1,1) = str2double(C{1,3});
    for jj = 2:NELEM
        while true
            current_line = current_line + 1;
            if ~(strcmp(lines{current_line}(1), '*') || strcmp(lines{current_line}(1), 'C'))
                break;
            end

        end
        C = strsplit(strtrim(lines{current_line,1}));
        IELEM(jj,1) = str2double(C{1,1});
    end
    DIV_GRP(ii,1).ID = ID;
    DIV_GRP(ii,1).NELEM = NELEM;
    DIV_GRP(ii,1).IELEM = IELEM;
end