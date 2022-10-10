function prop = readIWFM_StreamProp(filename, Nnodes, Nprop, Nskip)
    % read the entire file
    str = fileread(filename);
    % split it into lines
    lines = regexp(str, '\r\n|\r|\n', 'split')';
    
    prop = nan(Nnodes,Nprop);
    for ii = 1:Nnodes
        C = strsplit(strtrim(lines{Nskip+ii}));
        C = str2double(C);
        C(:,isnan(C)) = [];
        prop(ii,:) = C(1:Nprop);
    end
end