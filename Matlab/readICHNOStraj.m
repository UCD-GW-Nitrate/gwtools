function S = readICHNOStraj(filename)
cc = strsplit(filename, '.');
if strcmp(cc{end},'h5')
    info = h5info(filename);
    PVA = h5read(info.Filename,[info.Name,info.Datasets(3).Name])';
    ESID = h5read(info.Filename,[info.Name,info.Datasets(1).Name])';
    ER = h5read(info.Filename,[info.Name,info.Datasets(2).Name]);
    
    S(size(ESID,1),1).Eid = [];
    S(size(ESID,1),1).Sid = [];
    S(size(ESID,1),1).p = [];
    S(size(ESID,1),1).v = [];
    S(size(ESID,1),1).t = [];
    S(size(ESID,1),1).ER = [];
    
    for ii = 1:size(ESID,1)
        idx = find(PVA(:,1) == ESID(ii,1) & PVA(:,2) == ESID(ii,2));
        S(ii,1).Eid = ESID(ii,1);
        S(ii,1).Sid = ESID(ii,2);
        [~, I] = sort(PVA(idx,3));
        S(ii,1).p = PVA(idx(I),4:6);
        S(ii,1).v = PVA(idx(I),7:9);
        S(ii,1).t = PVA(idx(I),10);
        S(ii,1).ER = ER(ii);
    end
else
    str = fileread(filename);
    lines = regexp(str, '\r\n|\r|\n', 'split')';
    new_strmline = true;
    idx = 0;

    S(100000,1).Eid = [];
    S(100000,1).Sid = [];
    S(100000,1).p = [];
    S(100000,1).v = [];
    S(100000,1).t = [];
    S(100000,1).ER = [];

    for ii = 1:length(lines)
        if isempty(lines{ii,1})
            continue;
        end
        C = strsplit(lines{ii,1});
        pid = str2double(C{1});
        if pid == -9
            Eid = str2double(C{2});
            Sid = str2double(C{3});
            %[Eid Sid]
            S(idx,1).Eid = str2double(C{2});
            S(idx,1).Sid = str2double(C{3});
            S(idx,1).p = p(1:cnt_pid-1,:);
            S(idx,1).v = v(1:cnt_pid-1,:);
            S(idx,1).t = t(1:cnt_pid-1,:);
            S(idx,1).ER = C{4};
            new_strmline = true;
        else
            if new_strmline
                cnt_pid = 1;
                p = nan(100000,3);
                v = nan(100000,3);
                t = nan(100000,1);
                idx = idx + 1;
                new_strmline = false;
            end
            tmp = cellfun(@str2double,C);
            p(cnt_pid,:) = tmp(4:6);
            v(cnt_pid,:) = tmp(7:9);
            t(cnt_pid,:) = tmp(10);
            cnt_pid = cnt_pid + 1;
        end

    end
    S(idx+1:end,:) = [];
end
for ii = 1:length(S)
    S(ii,1).Xend = S(ii,1).p(end,1); 
    S(ii,1).Yend = S(ii,1).p(end,2);
    S(ii,1).Time = abs(S(ii,1).t(end));
    S(ii,1).Len = sum(sqrt(sum(diff(S(ii,1).p).^2,2)));
end

