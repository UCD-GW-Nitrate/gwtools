function mesh = readIWFM_Elements(filename, Nel, Nheader)
%readIWFM_Elements(filename, Nel, Nheader)
% Reads the elements of the IWFM mesh
% Inputs:
%   Nel is the number of elements
%   Nheader are the number of lines to skip before start reading the mesh
% indices. The Nheader should skip the region names.
fid = fopen(filename,'r');
temp = textscan(fid, '%f %f %f %f %f %f', Nel, 'HeaderLines', Nheader);
fclose(fid);
mesh = cell2mat(temp);
end

