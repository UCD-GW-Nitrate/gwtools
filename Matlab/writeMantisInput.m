function writeMantisInput(opt)
%writeMantisInput writes the options to a file to be used by the testclient
% program

fid = fopen(opt.infile,'w');
if ~isempty(opt.descr)
    for ii = 1:length(opt.descr)
        fprintf(fid, '# %s\n',opt.descr{ii})
    end
end

fprintf(fid, 'endSimYear %d\n', opt.endSimYear);
fprintf(fid, 'startRed %d\n', opt.startRed);
fprintf(fid, 'endRed %d\n', opt.endRed);
fprintf(fid, 'flowScen %s\n', opt.flowScen);
fprintf(fid, 'loadScen %s\n', opt.loadScen);
fprintf(fid, 'unsatScen %s\n', opt.unsatScen);
fprintf(fid, 'unsatWC %.2f\n', opt.unsatWC);
fprintf(fid, 'bMap %s\n', opt.bMap);

fprintf(fid, 'Nregions %d', length(opt.Regions));
for ii = 1:length(opt.Regions)
    fprintf(fid, ' %s',opt.Regions{ii});
end
fprintf(fid, '\n');

if ~isempty(opt.RadSelect)
    fprintf(fid, 'RadSelect %.2f %.2f %.2f\n', opt.RadSelect);
end
if ~isempty(opt.RectSelect)
    fprintf(fid, 'RectSelect %.2f %.2f %.2f %.2f\n', opt.RectSelect);
end
if ~isempty(opt.DepthRange)
    fprintf(fid, 'DepthRange %.2f %.2f\n', opt.DepthRange);
end
if ~isempty(opt.ScreenLenRange)
    fprintf(fid, 'ScreenLenRange %.2f %.2f\n', opt.ScreenLenRange);
end
if ~isempty(opt.SourceArea)
    fprintf(fid, 'SourceArea %d %d %d %.2f\n', opt.SourceArea);
end
if isempty(opt.Crops)
    opt.Crops = [-9 1];
end
fprintf(fid, 'Ncrops %d\n', size(opt.Crops,1));
fprintf(fid, '%d %.3f\n', opt.Crops');

if ~isempty(opt.PixelRadius)
    fprintf(fid, 'PixelRadius %d\n', opt.PixelRadius);
end

if ~isempty(opt.DebugID)
    fprintf(fid, 'DebugID %s\n', opt.DebugID);
end

fclose(fid);
end

