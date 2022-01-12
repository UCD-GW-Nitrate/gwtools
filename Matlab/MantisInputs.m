function opt = MantisInputs()
%MantisOptions Returns a structure with the simulation inputs
%   Detailed explanation goes here
opt.infile = 'incomingMSG.dat';
opt.outfile = 'testClientResults.dat';
opt.descr = {'This is a description of the simulation input','It will be ignored'};
opt.endSimYear = 2100;
opt.startRed = 2020;
opt.endRed = 2030;
opt.flowScen = 'C2VSIM_II_01';
opt.unsatScen = 'C2VSIM_SPRING_2000';
opt.wellType = 'VI';
opt.unsatWC = 0.01;
opt.rchScen = opt.flowScen;
opt.bMap = 'Subregions';
opt.Regions = {'Subregion1','Subregion2'};
opt.RadSelect = []; % X Y radius
opt.RectSelect = []; %Xmin Ymin Xmax Ymax
opt.DepthRange = []; %min max
opt.ScreenLenRange = []; %min max
opt.loadScen = 'GNLM';
opt.loadSubScen = '';
opt.modifierName = '';
opt.modifierType = '';
opt.isLoadConc = 1;
opt.Crops = [-9 1];
opt.SourceArea = []; %Npixels minPixels maxPixels percPixels
opt.PixelRadius = [];
opt.DebugID = [];
end

