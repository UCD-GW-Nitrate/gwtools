function xyz = readNPSAT_xyz(prefix, Nproc, Zeropad)
    xyz = [];
    frmt = ['%0' num2str(Zeropad) 'd'];
    for ii = 1:Nproc
        fid = fopen([prefix num2str(ii-1,frmt) '.xyz'],'r');
        Nc = cell2mat(textscan(fid, '%f',1));
        if Nc > 0
            C = cell2mat(textscan(fid, '%f %f %f %f',Nc));
            xyz = [xyz;C];
        end
        fclose(fid);
    end
end