function BUD = readIWFM_Zbudget(filename, Ntimes, Nzones)
if isempty(Ntimes)
    Ntimes = 1056;
end
if isempty(Nzones)
    Nzones = 1392;
end

fid = fopen(filename, 'r');
cnt_zones = 1;
for izone = 1:Nzones
    for ii = 1:2
        temp = fgetl(fid);
    end
    C = strsplit(temp);
    ZoneID = str2double(C{end});
    display(ZoneID);
    for ii = 1:4;temp = fgetl(fid);end
    labels = cellfun(@strtrim, strsplit(temp, '          '), 'UniformOutput', false);
    for ii = 1:2;temp = fgetl(fid);end
    DATA = nan(Ntimes, (length(labels)-2)*2+1);
    for ii = 1:Ntimes
        temp = fgetl(fid);
        C = strsplit(temp,' ');
        c = textscan(C{1,1},'%f/%f/%f/_%s');
        Year(ii,1) = c{1,3};
        Mon(ii,1) = c{1,1};
        C(:,1) =  [];
        DATA(ii,:) = cellfun(@str2double,C);
    end
    BUD(cnt_zones,1).ZoneID = ZoneID;
    BUD(cnt_zones,1).Headres = labels';
    BUD(cnt_zones,1).Data = DATA;
    cnt_zones = cnt_zones + 1;
    for ii = 1:2;temp = fgetl(fid);end
end

fclose(fid);