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
    error('Not Implemented yet')
    prefix = filename(1:dot_pos(end));
    

else
    error([ext ' is unknown extension'])
end


end