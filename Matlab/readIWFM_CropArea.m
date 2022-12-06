function CropAreas = readIWFM_CropArea(filename, Ntimes, Nelem, Ncrops, Nskip)
% CropAreas = readIWFM_CropArea(filename, Ntimes, Nelem, Nskip)
% Reads the crop area per element

% read the entire file
str = fileread(filename);
% split it into lines
lines = regexp(str, '\r\n|\r|\n', 'split')';

current_line = Nskip + 1;
CropAreas.Data{Ntimes,1} = zeros(Nelem, Ncrops+1);
CropAreas.YMD = zeros(Ntimes, 3);
for ii = 1:Ntimes
    C = strsplit(lines{current_line,1});
    cctime = textscan(C{1,1},'%f/%f/%f/_%s');
    CropAreas.YMD(ii,:) = [cctime{3} cctime{1} cctime{2}];
    C(:,1) = [];
    CropAreas.Data{ii,1} = zeros(Nelem, Ncrops+1);
    CropAreas.Data{ii,1}(1,:) = cellfun(@str2double,C);
    current_line = current_line + 1;
    for j = 2:Nelem
         C = strsplit(lines{current_line,1});
         CC = cellfun(@str2double,C);
         CC(:,isnan(CC)) = [];
         CropAreas.Data{ii,1}(j,:) = cellfun(@str2double,C);
         current_line = current_line + 1;
    end
end

