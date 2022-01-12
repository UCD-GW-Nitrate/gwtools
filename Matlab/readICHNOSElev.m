function Elev = readICHNOSElev(filename)
% Elev = readICHNOSElev(filename)
%
% Reads the files that describe the top and bottom boundary of the ICHNOS
% code

str = fileread(filename);
lines = regexp(str, '\r\n|\r|\n', 'split')';

idx = 1;
% Read the interpolation type
C = textscan(lines{idx,1},'%s');
idx = idx + 1;
Elev.Type = C{1}{1};

if Elev.Type == 'CLOUD'
    % Read the parameters
    C = textscan(lines{idx,1},'%f');
    idx = idx + 1;
    Elev.R = C{1}(1);
    Elev.P = C{1}(2);
    pp = [];
    while 1
        C = textscan(lines{idx,1},'%f');
        idx = idx + 1;
        pp = [pp;C{1}'];
        if idx > length(lines)
        break;
        end
        if isempty(lines{idx,1})
           break; 
        end
    end
    Elev.Data = pp;
elseif Elev.Type == 'MESH2D'
    warning('Reading MESH2D is not implemented yet');
end