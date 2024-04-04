function VIRTUAL_WELLS_CLASS = MantisAssignWells2Regions(VIRTUAL_WELLS,mantis_path)
% VIRTUAL_WELLS_CLASS = MantisAssignWells2Regions(VIRTUAL_WELLS,mantis_path)

% Read the background maps
basins = shaperead(fullfile(mantis_path,'gis_data','BackgroundMaps','Basins',"C2VsimBasins_3310"));
counties = shaperead(fullfile(mantis_path,'gis_data','BackgroundMaps','Counties','counties_simple_CVclip'));
b118 = shaperead(fullfile(mantis_path,'gis_data','BackgroundMaps','B118','B118_CV_TAproj'));
townships = shaperead(fullfile(mantis_path,'gis_data','BackgroundMaps','Townships','CVHM_Townships_3310_simplified'));
subreg = shaperead(fullfile(mantis_path,'gis_data','BackgroundMaps','Subregions','C2Vsim_Subregions_3310'));

% Find in which areas the wells belong to
VIRTUAL_WELLS_CLASS = table('Size',[size(VIRTUAL_WELLS,1),5], 'VariableTypes',{'categorical','categorical','categorical','categorical','categorical'},...
    'VariableNames',{'Basins', 'Counties', 'B118', 'Townships', 'Subregions'});

% Basins
basinNames = categorical({'SanJoaquinValley', 'SacramentoValley',  'TulareLakeBasin'});
for ii = 1:length(basins)
    basins(ii,1).pp = polyshape(basins(ii,1).X, basins(ii,1).Y);
    %idx = find(my_inpolygon(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y, basins(ii,1).X, basins(ii,1).Y));
    idx = find(basins(ii,1).pp.isinterior(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y));
    VIRTUAL_WELLS_CLASS.Basins(idx) = basinNames(ii);
end

% For the wells that are not defined assign the category of the closest defined well
id_undef = find(isundefined(VIRTUAL_WELLS_CLASS.Basins));
id_def = find(~isundefined(VIRTUAL_WELLS_CLASS.Basins));
for ii = 1:length(id_undef)
    dst = sqrt((VIRTUAL_WELLS.X(id_undef(ii)) - VIRTUAL_WELLS.X(id_def)).^2 + (VIRTUAL_WELLS.Y(id_undef(ii)) - VIRTUAL_WELLS.Y(id_def)).^2);
    [c, d] = min(dst);
    VIRTUAL_WELLS_CLASS.Basins(id_undef(ii)) = VIRTUAL_WELLS_CLASS.Basins(id_def(d));
end

% Counties
for ii = 1:length(counties)
    classname = replace(counties(ii,1).name, " ", "");
    counties(ii,1).pp = polyshape(counties(ii,1).X, counties(ii,1).Y);
    %id = find(my_inpolygon(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y, counties(ii,1).X, counties(ii,1).Y));
    id = find(counties(ii,1).pp.isinterior(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y));
    VIRTUAL_WELLS_CLASS.Counties(id) = classname;
end

id_undef = find(isundefined(VIRTUAL_WELLS_CLASS.Counties));
id_def = find(~isundefined(VIRTUAL_WELLS_CLASS.Counties));
for ii = 1:length(id_undef)
    dst = sqrt((VIRTUAL_WELLS.X(id_undef(ii)) - VIRTUAL_WELLS.X(id_def)).^2 + (VIRTUAL_WELLS.Y(id_undef(ii)) - VIRTUAL_WELLS.Y(id_def)).^2);
    [c, d] = min(dst);
    VIRTUAL_WELLS_CLASS.Counties(id_undef(ii)) = VIRTUAL_WELLS_CLASS.Counties(id_def(d));
end

% B118
for ii = 1:length(b118)
    classname = replace(b118(ii,1). Basin_Subb,{' ','-','.'}, '_');
    b118(ii,1).pp = polyshape(b118(ii,1).X, b118(ii,1).Y);
    %id = find(my_inpolygon(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y, b118(ii,1).X, b118(ii,1).Y));
    id = find(b118(ii,1).pp.isinterior(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y));
    VIRTUAL_WELLS_CLASS.B118(id) = classname;
end

% For the wells that are not defined assign the category of the closest defined well
id_undef = find(isundefined(VIRTUAL_WELLS_CLASS.B118));
id_def = find(~isundefined(VIRTUAL_WELLS_CLASS.B118));
for ii = 1:length(id_undef)
    dst = sqrt((VIRTUAL_WELLS.X(id_undef(ii)) - VIRTUAL_WELLS.X(id_def)).^2 + (VIRTUAL_WELLS.Y(id_undef(ii)) - VIRTUAL_WELLS.Y(id_def)).^2);
    [c, d] = min(dst);
    VIRTUAL_WELLS_CLASS.B118(id_undef(ii)) = VIRTUAL_WELLS_CLASS.B118(id_def(d));
end

% Townships
for ii = 1:length(townships)
    townships(ii,1).pp = polyshape(townships(ii,1).X, townships(ii,1).Y);
    %id = find(my_inpolygon(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y, townships(ii,1).X, townships(ii,1).Y));
    id = find(townships(ii,1).pp.isinterior(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y));
    VIRTUAL_WELLS_CLASS.Townships(id) = townships(ii,1).CO_MTR;
end

id_undef = find(isundefined(VIRTUAL_WELLS_CLASS.Townships));
id_def = find(~isundefined(VIRTUAL_WELLS_CLASS.Townships));
for ii = 1:length(id_undef)
    dst = sqrt((VIRTUAL_WELLS.X(id_undef(ii)) - VIRTUAL_WELLS.X(id_def)).^2 + (VIRTUAL_WELLS.Y(id_undef(ii)) - VIRTUAL_WELLS.Y(id_def)).^2);
    [c, d] = min(dst);
    VIRTUAL_WELLS_CLASS.Townships(id_undef(ii)) = VIRTUAL_WELLS_CLASS.Townships(id_def(d));
end

% Sub-basins
for ii = 1:length(subreg)
    subreg(ii,1).pp = polyshape(subreg(ii,1).X, subreg(ii,1).Y);
    %id = find(my_inpolygon(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y, subreg(ii,1).X, subreg(ii,1).Y));
    id = find(subreg(ii,1).pp.isinterior(VIRTUAL_WELLS.X, VIRTUAL_WELLS.Y));
    VIRTUAL_WELLS_CLASS.Subregions(id) = categorical({['Subregion' num2str(subreg(ii,1).IRGE)]});
end

% For the wells that are not defined assign the category of the closest defined well
id_undef = find(isundefined(VIRTUAL_WELLS_CLASS.Subregions));
id_def = find(~isundefined(VIRTUAL_WELLS_CLASS.Subregions));
for ii = 1:length(id_undef)
    dst = sqrt((VIRTUAL_WELLS.X(id_undef(ii)) - VIRTUAL_WELLS.X(id_def)).^2 + (VIRTUAL_WELLS.Y(id_undef(ii)) - VIRTUAL_WELLS.Y(id_def)).^2);
    [c, d] = min(dst);
    VIRTUAL_WELLS_CLASS.Subregions(id_undef(ii)) = VIRTUAL_WELLS_CLASS.Subregions(id_def(d));
end
end