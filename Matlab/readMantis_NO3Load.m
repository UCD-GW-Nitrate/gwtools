function LD = readMantis_NO3Load(filename)
%LD = readMantis_NO3Load(filename) Reads the NO3 loading files
%
% If the data are printed in h5 file then filename must include the
% extension h5. If the loading is written in the ASCII format which
% consists of 2 files use either the idxlu or the nload

dot_pos = find(filename == '.');
ext = filename((dot_pos(end)+1):end);

if strcmp(ext,'h5')
    LDinfo = h5info(filename);
    for ii = 1:length(LDinfo.Datasets)
        LD.(LDinfo.Datasets(ii).Name) = h5read(LDinfo.Filename,[LDinfo.Name,LDinfo.Datasets(ii).Name]);
    end
elseif strcmp(ext,'idxlu') || strcmp(ext,'nload')
    prefix = filename(1:dot_pos(end));
    % Read the idxlu
    fid = fopen([prefix 'idxlu'],'r');
    C = textscan(fid,'%d %d %d %d',1);
    Ncells = C{1,1};
    LD.LUNgrid = zeros(1,6);
    LD.LUNgrid(1) = double(C{1,2});
    LD.LUNgrid(2) = double(C{1,3});
    LD.LUNgrid(3) = double(C{1,4});
    C = textscan(fid, ['%d ' repmat('%d', 1,LD.LUNgrid(3))], Ncells);
    LD.Nidx = C{1,1};
    LD.LU = nan(Ncells, LD.LUNgrid(3));
    for ii = 1:LD.LUNgrid(3)
        LD.LU(:,ii) = C{1,ii+1};
    end
    fclose(fid);

    %Read the nload
    fid = fopen([prefix 'nload'],'r');
    C = textscan(fid,'%d %d %d %d',1);
    Nts = C{1,1};
    LD.LUNgrid(4) = double(C{1,2});
    LD.LUNgrid(5) = double(C{1,3});
    LD.LUNgrid(6) = double(C{1,4});
    C = textscan(fid, ['%f ' repmat('%f', 1,LD.LUNgrid(6))], Nts);
    
    LD.Nload = nan(Nts, LD.LUNgrid(6));
    for ii = 1:LD.LUNgrid(6)
        LD.Nload(:,ii) = C{1,ii};
    end
    fclose(fid);
else
    error([ext ' is unknown extension'])
end


end