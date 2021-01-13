function h = readIWFM_headalloutput(varargin)

if nargin == 1
    filename = varargin{1};
    printout = true;
end

if nargin == 2
    filename = varargin{1};
    printout = varargin{2};
end

fid = fopen(filename, 'r');
cnt_per = 0;
cnt_lay = 1;
h{1,2} = [];

    
while 1 
   try 
      temp = fgetl(fid); 
      if isempty(temp)
          continue
      end
      if strcmp(temp(1),'*')
          continue
      end
      C = strsplit(temp,' ');
      for ii = 1:length(C)
          if isempty(C{1,ii})
              continue
          end
          c = textscan(C{1,ii},'%f/%f/%f/_%s');
          if isempty(c{1,2})
              %Then its a head value
              h{cnt_per,2}(cnt_nodes, cnt_lay) = c{1,1};
              cnt_nodes = cnt_nodes + 1;
              

          else
              %its time stamp
              cnt_per = cnt_per + 1;
              h{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
              if printout
                display(h{cnt_per,1})
              end
              cnt_nodes = 1;
              cnt_lay = 1;
          end
      end
      cnt_lay = cnt_lay + 1;
      cnt_nodes = 1;
   catch
       break;
   end
end
fclose(fid);

