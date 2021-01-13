function VelOut = readIWFM_Velocity(filename, Nel, Nlay)
%readC2Vsim_Velocity Read the velocity field that is specified in the 
% VELOUTFL parameter
% Nel is the number of elements which is 32537 for the C2VsimFGv1 version
% Nlay is the number of layers which is 4

%   Detailed explanation goes here

str = fileread(filename);
lines = regexp(str, '\r\n|\r|\n', 'split')';
section = 1;
ND = nan(Nel,2);
idx = 0;
idx_step = 1;
YMD = nan(2000,3);
VX = nan(Nel, 2000, Nlay);
VY = nan(Nel, 2000, Nlay);
VZ = nan(Nel, 2000, Nlay);

while idx < length(lines)
    idx = idx+1;
    if isempty(lines{idx,1})
        continue;
    end
    if strcmp(lines{idx,1}(1),'*')
        continue;
    end
    if section == 1
        for j = 1:Nel
            cc = textscan(lines{idx,1},'%f');
            idx = idx + 1;
            ND(cc{1}(1),:) = [cc{1}(2) cc{1}(3)];
        end
        section = section + 1;
        continue;
    end
    if section == 2
        C = strsplit(lines{idx,1},' ');
        display(C{1,1});
        c = textscan(C{1,1},'%f/%f/%f/_%s');
        YMD(idx_step,:) = [c{1,3} c{1,1} c{1,2}];
        C(:,1) =  [];
        tmp = cellfun(@str2double,C)';
        tmp1 = reshape(tmp(2:end), length(tmp(2:end))/Nlay, Nlay);
        VX(tmp(1),idx_step,:) = tmp1(1,:);
        VY(tmp(1),idx_step,:) = tmp1(2,:);
        VZ(tmp(1),idx_step,:) = tmp1(3,:);
        idx = idx + 1;
        for j = 2:Nel
            cc = cell2mat(textscan(lines{idx,1},'%f'));
            tmp = reshape(cc(2:end), length(cc(2:end))/Nlay, Nlay);
            VX(cc(1),idx_step,:) = tmp(1,:);
            VY(cc(1),idx_step,:) = tmp(2,:);
            VZ(cc(1),idx_step,:) = tmp(3,:);
            if j ~= Nel
                idx = idx + 1;
            end
        end
        idx_step = idx_step + 1;
        
    end

end
YMD(idx_step:end,:) = [];
VX(:,idx_step:end,:) = [];
VY(:,idx_step:end,:) = [];
VZ(:,idx_step:end,:) = [];
VelOut.YMD = YMD;
VelOut.ND = ND;
VelOut.VX = VX;
VelOut.VY = VY;
VelOut.VZ = VZ;
end

