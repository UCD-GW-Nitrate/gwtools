function h = readIWFM_StreamStage(filename, Ntimes)
% h = readIWFM_StreamStage(filename, Ntimes)
%
% Reads the Streamstage filethat is specified in the 
% STHYDOUTFL parameter under the stream main file.
%
% Keep in mind that depending on the flag IHSQR the output of this file can
% be stream flows or stream stage.
%
% Ntimes is the number of time steps in the simulation. By default this is
% 504

str = fileread(filename);
lines = regexp(str, '\r\n|\r|\n', 'split')';
for ii = 1:length(lines)
    if ~(strcmp(lines{ii}(1), '*') || strcmp(lines{ii}(1), 'C'))
        break;
    end
end
iheader = ii - 3;

% Read Hydrograph IDs
C = strsplit(lines{iheader,1});
C(:,1:3) = [];
HydId = cellfun(@str2double,C)';
C = strsplit(lines{iheader+1,1});
C(:,1:2) = [];
IRV = cellfun(@str2double,C)';

iskip = iheader + 2;
Hs = nan(length(IRV), Ntimes);
for ii = 1:Ntimes
    C = strsplit(lines{iskip + ii,1});
    cctime = textscan(C{1,1},'%f/%f/%f/_%s');
    YMD(ii,:) = [cctime{3} cctime{1} cctime{2}];
    C(:,1) = [];
    Hs(:,ii) = cellfun(@str2double,C)';
end
h.HydId = HydId;
h.IRV = IRV;
h.YMD = YMD;
h.Hs = Hs;
end

