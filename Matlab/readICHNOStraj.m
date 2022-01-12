function S = readICHNOStraj(filename)

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