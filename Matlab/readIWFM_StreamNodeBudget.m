function BUD = readIWFM_StreamNodeBudget(filename, nTimeSteps)
    strm_node_bud = h5info(filename);
    Stream_col_names = strm_node_bud.Groups.Attributes(11).Value(2:end);
    for j = 1:length(Stream_col_names)
        DATA{j,1} = zeros(length(strm_node_bud.Datasets), nTimeSteps);
    end
    

    for ii = 1:length(strm_node_bud.Datasets)
        id = textscan(strm_node_bud.Datasets(ii).Name,'NODE%s');
        id = str2double(id{1,1}{1,1});
        ttt = h5read(strm_node_bud.Filename,[strm_node_bud.Name strm_node_bud.Datasets(ii).Name]);
        for j = 1:length(Stream_col_names)
            DATA{j,1}(id,:) = ttt(j, :);
        end
    end
    
    BUD.ColNames = Stream_col_names;
    BUD.Data = DATA;
end