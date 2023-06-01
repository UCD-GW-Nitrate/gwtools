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

% -----------------
% Mantis client options
% -----------------
opt.infile = 'incomingMSG.dat';
opt.outfile = 'testClientResults.dat';
opt.descr = {'This is a description of the simulation input','It will be ignored'};

% -----------------
% GUI Scenario Options
% -----------------
opt.modelArea = 'CentralValley';
opt.bMap = 'Townships';
opt.Regions = {'39M02S10E','39M02S09E'};
opt.flowScen = 'C2VSIM_III_03';
opt.wellType = 'VI';
opt.unsatScen = 'C2VSIM_SPRING_2000';
opt.unsatWC = 0.01;
opt.endSimYear = 2100;

% -----------------
% loading options
% -----------------
opt.loadScen = 'SWAT1';
opt.startRed = 2020;
opt.endRed = 2030;
opt.Crops =[-9 0.85];

% -----------------
% Options that take default values
% -----------------
opt.minRch = 10;
opt.maxConc = 250;
opt.constRed = 1.0;
opt.LoadTransitionName = 'GNLM';
opt.LoadTransitionStart = 2005;
opt.LoadTransitionEnd = 2015;
opt.startSimYear = 1925;

% -----------------
% Filter Options
% -----------------
opt.RadSelect = [];
opt.RectSelect = [];
opt.DepthRange = [];
opt.ScreenLenRange = [];
opt.SourceArea = [];

% -----------------
% Raster loading Options
% -----------------
opt.loadSubScen = '';
opt.modifierName = '';
opt.modifierType = '';
opt.modifierUnit = '';

% -----------------
% Misc Options
% -----------------
opt.getids = [];
opt.DebugID = [];
opt.printLF = [];
opt.printURF = [];
opt.printBTC = [];
opt.printWellBTC = [];
opt.PixelRadius = [];

end

