function h = readIWFM_headalloutput(varargin)
% readIWFM_headalloutput reads the output of the outputfile specified at
% GWALLOU in the main groundwater input file.
% This function can be run with two modes 
%
% One is relatively slow but works without specifying additional
% inputs such as number of nodes layers etc. It should identify everything
% from the head output file.
%
% USAGE
% h = readIWFM_headalloutput(filename)
% h = readIWFM_headalloutput(filename, 0 or 1)
%
% Inputs: 
% filename
%
% The second mode is a bit faster but requires the following inputs
%
% USAGE
% h = readIWFM_headalloutput(filename, Nnodes, Nlay, Ntimes, 0 or 1)
%
% Inputs: 
% filename
% Nnodes: Number of nodes
% Nlay: Number of layers
% Ntimes: Number of times
%
% In both modes you can specify to print or not the reading progress by
% setting
% printout : 0 not print, 1 print 


if nargin == 1
    filename = varargin{1};
    h = headReadV1(filename, 0);
elseif nargin == 2
    filename = varargin{1};
    printout = varargin{2};
    h = headReadV1(filename, printout);
elseif nargin == 5
    filename = varargin{1};
    Nnodes = varargin{2};
    Nlay = varargin{3};
    Ntimes = varargin{4};
    printout = varargin{5};
    h = headReadV2(filename, Nnodes, Nlay, Ntimes, printout);
else
    error('Wrong number of arguments. They must be 1,2 or 5')
end

end

function h = headReadV2(filename, Nnodes, Nlay, Ntimes, printout)
    % read the entire file
    str = fileread(filename);
    % split it into lines
    lines = regexp(str, '\r\n|\r|\n', 'split')';
    
    cnt_per = 0;
    cnt_lay = 1;
    h{Ntimes,2} = [];
    ilay = 1;
    
    for iline = 1:length(lines)
    
        if isempty(lines{iline,1})
            continue
        end
    
        if strcmp(lines{iline,1}(1),'*')
            continue
        end
        C = strsplit(lines{iline,1},' ');
        if ilay == 1
            c = textscan(C{1,1},'%f/%f/%f/_%s');
            cnt_per = cnt_per + 1;
            h{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
            h{cnt_per,2} = nan(Nnodes,Nlay);
            if printout
                display(h{cnt_per,1})
    
            C(:,1) = [];
            end
        end
        C = str2double(C);
        C(:,isnan(C)) = [];
        h{cnt_per,2}(:, ilay) = C';
        ilay = ilay + 1;
        if ilay > Nlay
            ilay = 1;
        end
    end
end


function h = headReadV1(filename, printout)
    fid = fopen(filename, 'r');
    cnt_per = 0;
    cnt_lay = 1;
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
end
