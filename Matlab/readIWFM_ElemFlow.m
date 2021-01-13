function F = readIWFM_ElemFlow(filename)
fid = fopen(filename, 'r');

cnt_lines = 0;
while 1 
   temp = fgetl(fid);
   if isempty(temp)
       continue
   end
   if temp ==-1
       break;
   end
   cnt_lines = cnt_lines + 1;
   if cnt_lines == 5
      % Read the layers
      C = strsplit(temp);
      C(:,1:2) =  [];
      LAY = cellfun(@str2double,C)';
      LAY(isnan(LAY),:) = [];
      temp = fgetl(fid);
      C = strsplit(temp);
      C(:,1:2) =  [];
      ELEM = nan(length(LAY),2);
      for ii = 1:length(C)
          if isempty(C{1,ii})
              continue;
          end
          ELEM(ii,:) = double(cell2mat(textscan(C{1,ii},'E%d-E%d')));
      end
      temp = fgetl(fid);
      DATA = nan(length(LAY),2000);
      Year = nan(2000,1);
      Mon = nan(2000,1);
      cnt_tm = 1;
   end
   
   if cnt_lines >= 6
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
F.Lay = LAY;
F.Year = Year(1:cnt_tm-1,:);
F.Mon = Mon(1:cnt_tm-1,:);
F.Elem = ELEM;
F.Flow = DATA(:,1:cnt_tm-1);
fclose(fid);

