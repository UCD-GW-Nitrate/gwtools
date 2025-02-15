function writeMantisInput(opt)
%writeMantisInput writes the options to a file to be used by the testclient
% program

fid = fopen(opt.infile,'w');
if ~isempty(opt.descr)
    for ii = 1:length(opt.descr)
        fprintf(fid, '# %s\n',opt.descr{ii});
    end
end

% -----------------
% GUI Scenario Options
% -----------------
fprintf(fid, 'modelArea %s\n', opt.modelArea);
fprintf(fid, 'bMap %s\n', opt.bMap);
fprintf(fid, 'Nregions %d', length(opt.Regions));
for ii = 1:length(opt.Regions)
    fprintf(fid, ' %s',opt.Regions{ii});
end
fprintf(fid, '\n');
fprintf(fid, 'flowScen %s\n', opt.flowScen);
fprintf(fid, 'wellType %s\n', opt.wellType);
fprintf(fid, 'unsatScen %s\n', opt.unsatScen);
fprintf(fid, 'unsatWC %.2f\n', opt.unsatWC);
fprintf(fid, 'endSimYear %d\n', opt.endSimYear);

% -----------------
% loading options
% -----------------
fprintf(fid, 'loadScen %s\n', opt.loadScen);
fprintf(fid, 'startRed %d\n', opt.startRed);
fprintf(fid, 'endRed %d\n', opt.endRed);
if isempty(opt.Crops)
    opt.Crops = [-9 1];
end
fprintf(fid, 'Ncrops %d\n', size(opt.Crops,1));
fprintf(fid, '%d %.3f\n', opt.Crops');

% -----------------
% Options that take default values
% -----------------
if ~isempty(opt.minRch)
    fprintf(fid, 'minRch %.2f\n', opt.minRch);
end
if ~isempty(opt.rchMap)
    fprintf(fid, 'rchMap %s\n', opt.rchMap);
end
if ~isempty(opt.maxConc)
    fprintf(fid, 'maxConc %.2f\n', opt.maxConc);
end
if ~isempty(opt.maxAge)
    fprintf(fid, 'maxAge %.2f\n', opt.maxAge);
end
if ~isempty(opt.constRed)
    fprintf(fid, 'constRed %.3f\n', opt.constRed);
end
if ~isempty(opt.LoadTransitionName)
    if isempty(opt.LoadTransitionStart) || isempty(opt.LoadTransitionEnd)
        opt.LoadTransitionStart = 2005;
        opt.LoadTransitionEnd = 2015;
    end
    fprintf(fid, 'loadTrans %s %d %d\n', opt.LoadTransitionName,...
        opt.LoadTransitionStart, opt.LoadTransitionEnd);
end
if ~isempty(opt.startSimYear)
    fprintf(fid, 'startSimYear %d\n', opt.startSimYear);
end
if ~isempty(opt.por)
    fprintf(fid, 'por %d\n', opt.por);
end
if ~isempty(opt.urfType)
    fprintf(fid, 'urfType %s\n', opt.urfType);
end

if ~isempty(opt.ADELR)
    fprintf(fid, 'ADELR %.5e %.5e\n', opt.ADELR);
end

if ~isempty(opt.initConcParam)
    fprintf(fid, 'initConcParam %.5e %.5e %5e %5e\n', opt.initConcParam);
end

% -----------------
% Filter Options
% -----------------
if ~isempty(opt.RadSelect)
    fprintf(fid, 'RadSelect %.2f %.2f %.2f\n', opt.RadSelect);
end
if ~isempty(opt.RectSelect)
    fprintf(fid, 'RectSelect %.2f %.2f %.2f %.2f\n', opt.RectSelect);
end
if ~isempty(opt.DepthRange)
    fprintf(fid, 'DepthRange %.2f %.2f\n', opt.DepthRange);
end
if ~isempty(opt.UnsatRange)
    fprintf(fid, 'UnsatRange %.2f %.2f\n', opt.UnsatRange);
end
if ~isempty(opt.Wt2tRange)
    fprintf(fid, 'Wt2tRange %.2f %.2f\n', opt.Wt2tRange);
end
if ~isempty(opt.ScreenLenRange)
    fprintf(fid, 'ScreenLenRange %.2f %.2f\n', opt.ScreenLenRange);
end
if ~isempty(opt.SourceArea)
    fprintf(fid, 'SourceArea %d %d %d %.2f\n', opt.SourceArea);
end

% -----------------
% Raster loading Options
% -----------------
if ~isempty(opt.loadSubScen)
    fprintf(fid, 'loadSubScen %s\n', opt.loadSubScen);
end
if ~isempty(opt.modifierName)
    fprintf(fid, 'modifierName %s\n', opt.modifierName);
end
if ~isempty(opt.modifierType)
    fprintf(fid, 'modifierType %s\n', opt.modifierType);
end
if ~isempty(opt.modifierUnit)
    fprintf(fid, 'modifierUnit %s\n', opt.modifierUnit);
end

% -----------------
% Misc Options
% -----------------
if ~isempty(opt.getids)
    if opt.getids ~= 0
        fprintf(fid, 'getids 1\n');
        if ~isempty(opt.getOtherInfo)
            if opt.getOtherInfo ~= 0
                fprintf(fid, 'getotherinfo 1\n');
            end
        end
    end
end
if ~isempty(opt.DebugID)
    fprintf(fid, 'DebugID %s\n', opt.DebugID);
    if ~isempty(opt.printLF)
        fprintf(fid, 'printLF %d\n', opt.printLF);
    end
    if ~isempty(opt.printURF)
        fprintf(fid, 'printURF %d\n', opt.printURF);
    end
    if ~isempty(opt.printBTC)
        fprintf(fid, 'printBTC %d\n', opt.printBTC);
    end
    if ~isempty(opt.printWellBTC)
        fprintf(fid, 'printWellBTC %d\n', opt.printWellBTC);
    end
end
if ~isempty(opt.PixelRadius)
    fprintf(fid, 'PixelRadius %d\n', opt.PixelRadius);
end

if ~isempty(opt.maxSourceCells)
    fprintf(fid, 'maxSourceCells %d\n', opt.maxSourceCells);
end


fclose(fid);
end

