function SmallWS = readIWFM_SmallWatersheds(filename, NSW, headerlines)
%SmallWS = readIWFM_SmallWatersheds(filename, NSW, headerlines)
% 
% Read the geometry info from the small watershed file. 
%
%   filename is the name of the file
%   NSW is the number of the watersheds
%   headerlines is the number of lines before the watershed listing

SmallWS(NSW,1).ID = [];
SmallWS(NSW,1).Areas = [];
SmallWS(NSW,1).IWBTS = [];
SmallWS(NSW,1).IWB = [];

fid = fopen(filename,'r');



for ii = 1:headerlines
    tmp = fgetl(fid);
end

for ii = 1:NSW
    temp = textscan(fid, '%f %f %f %f %f %f', 1);
    SmallWS(ii,1).ID = temp{1,1}(1);
    SmallWS(ii,1).Areas = temp{1,2}(1);
    SmallWS(ii,1).IWBTS = temp{1,3}(1);
    n = temp{1,4}(1);
    IWB = zeros(n,1);
    IWB(1,1) = temp{1,5}(1);
    if n > 1
        temp = textscan(fid, '%f %f', n-1);
        IWB(2:end,1) = temp{1,1};
    end
    SmallWS(ii,1).IWB = IWB;
end

fclose(fid);