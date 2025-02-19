function LD = readMantisLinearData(filename)
%readMantisLinearData read the data from the linear data format
C = strsplit(filename,'.');
if strcmp(C{end},'h5')
    LDinfo = h5info(filename);
    for ii = 1:length(LDinfo.Datasets)
        LD.(LDinfo.Datasets(ii).Name) = h5read(LDinfo.Filename,[LDinfo.Name,LDinfo.Datasets(ii).Name]);
    end
else
    fid = fopen(filename,'r');
    Ndata = cell2mat(textscan(fid,'%f',1));
    for ii = 1:Ndata
        tmp = textscan(fid,'%s',1);
        LD.Names{ii,1} = tmp{1};
    end
    fclose(fid);
    frmt = '%f';
    if Ndata > 1
        frmt = [frmt repmat(' %f',1,Ndata-1)];
    end
    fid = fopen(filename,'r');
    C = textscan(fid, frmt,'HeaderLines',Ndata+1);
    fclose(fid);
    LD.Data = cell2mat(C);
end
