function opt = MantisInputs(varargin)
% MantisOptions Returns a structure with the simulation inputs
%   It can be run without arguments such as 
%   scenario = MantisInputs()
%   In that case the scenario.client is empty and must be set before
%   running the MantisServere
%   or
%   scenario = MantisInputs(client)
%   where client is the path of the MantisClient. This sets the 
%   scenario.client = client


if ~isempty(varargin)
    opt.client = varargin{1};
else
    opt.client = '';
end

opt.infile = 'incomingMSG.dat';
opt.outfile = 'testClientResults.dat';
opt.descr = {'This is a description of the simulation input','It will be ignored'};
opt.endSimYear = 2100;
opt.startRed = 2020;
opt.constRed = [];
opt.endRed = 2030;
opt.flowScen = 'C2VSIM_II_03';
opt.unsatScen = 'C2VSIM_SPRING_2000';
opt.wellType = 'VI';
opt.unsatWC = 0.01;
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
opt.modifierUnit = '';
opt.Crops = [-9 1];
opt.minRch = 10;
opt.rchMap = [];
opt.maxConc = -1;
opt.loadTrans = '';
opt.loadTransYearStart = 2005;
opt.loadTransYearEnd = 2020;
opt.SourceArea = []; %Npixels minPixels maxPixels percPixels
opt.PixelRadius = [];
opt.DebugID = [];
opt.getids = [];
end

