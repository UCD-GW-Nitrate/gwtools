function out = readIWFM_VertFlow(filename)
fid = fopen(filename, 'r');
cnt_lines = 0;
cnt_tm = 1;
while 1 
    temp = fgetl(fid);
    if isempty(temp)
       continue
   end
   if temp ==-1
       break;
   end
   cnt_lines = cnt_lines + 1;
   if cnt_lines == 6
       C = strsplit(temp);
       C(:,1:2) = [];
       Region = cellfun(@str2double,C)';
       Region(isnan(Region),:) = [];
       for ii = 1:2;temp = fgetl(fid);end
       cnt_lines = cnt_lines + 1;
       DATA = nan(length(Region)*2,2000);
       Year = nan(2000,1);
       Mon = nan(2000,1);
   end
   if cnt_lines >= 8
       C = strsplit(temp,' ');
       c = textscan(C{1,1},'%f/%f/%f/_%s');
       display(C{1,1});
       Year(cnt_tm,1) = c{1,3};
       Mon(cnt_tm,1) = c{1,1};
       C(:,1) =  [];
       DATA(:,cnt_tm) = cellfun(@str2double,C)';
       cnt_tm = cnt_tm + 1;
   end
end
out.Regions = Region;
out.Year = Year(1:cnt_tm-1,:);
out.Mon = Mon(1:cnt_tm-1,:);
out.Flow = DATA(:,1:cnt_tm-1);
fclose(fid);