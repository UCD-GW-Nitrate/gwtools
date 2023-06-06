function out_raster = Equal_Areas(in_raster, Nunits, id, excludeid, method, opt)
    % out_raster = Equal_Areas(in_raster, Nunits, id, excludeid, method, opt)
    % This function divides a raster into equal areas
    %
    % INPUTS:
    % in_raster     : raster to be divided. Each cell id in the raster is
    %                 treated as a different region
    % Nunits        : The number of cells that each each equal area cell will have (Units are cells).
    % id            : Divides only the area associated with the particular id.
    %                 If empty then all areas are divided. 
    %
    % excludeid     : any area associated with ids in this list will not be divided
    %                   This is usefull to exclude no data cell.
    % method        : switches between the scott and metis method
    % opt           : options for the metis case. Consists of the following fields  
    %          name : name of the simulation.
    %    metis_path : Path of the metis program. It can be empty if
    %                 Default_metis_path variable is defined
    %
    % OUTPUTS
    % outraster     : The raster with the equal area cells.


switch method
    case 'scott'
        out_raster = Equal_Area_Scott(in_raster, Nunits, id, excludeid);
        
    case 'metis'
        All_ids = unique(in_raster);
        Xcoords = 1:size(in_raster,2);
        Ycoords = size(in_raster,1):-1:1;
        out_raster = nan(size(in_raster));
        for i = 1:length(All_ids)
            if ~isempty(find(excludeid == All_ids(i),1))
                continue
            end
            N_pixel_per_region = sum(sum(in_raster == All_ids(i)));
            N_eq_areas = round(N_pixel_per_region / Nunits);
            clear p TR;
            [ii, jj] = find(in_raster == All_ids(i));
            p(:,1) = Xcoords(jj);
            p(:,2) = Ycoords(ii);
            TR = delaunay( p(:, 1), p(:, 2) );
            %run metis
            if N_eq_areas < 2
                N_eq_areas = 2;
            end
            opt = check_metis_options(opt);
            [ElPart, NDPart] = decompose_2D_mesh(TR, N_eq_areas, opt);
            
            out_raster(sub2ind(size(out_raster), ii, jj)) = All_ids(i)*1000+NDPart+1;
        end
        
    otherwise
        warning('Uknown method. Use scott or metis');
        out_raster = nan;
        
end

end


function [ElPart NDPart]=decompose_2D_mesh(MSH_2d,Nsub,opt)
    % This function decompose a 2D mesh into Nsub subdomains
    %
    %Input:
    % p_2d :[np x 2] coordinates of mesh points
    %
    % MSH_2d : Mesh structure
    %
    % Nsub : either a number or a string. In case of string specifies the file
    % that decribes the subdomain division. Otherwise runs the Metis and create
    % the file.
    %
    % Nlay is the number of layers. The elevation has to be adjusted later

    % opt structure variable with fields 
    %    name : name of the simulation. In case number Nsub, this is the name 
    %    metis_path : Path of the metis program. It can be empty if
    %                 Default_metis_path variable is defined
    %

    %create inputfile for mpmetis
frmt='%d';
for ii = 2:size(MSH_2d,2)-1
    frmt = [frmt ' %d'];
end
frmt = [frmt ' %d\n'];
fid = fopen([opt.name '.mesh'],'w');
fprintf(fid,'%d\n',size(MSH_2d,1));
fprintf(fid,frmt,MSH_2d');
fclose(fid);

    %run metis
system([opt.metis_path ' ' opt.name '.mesh ' num2str(Nsub)])

    %read mesh partition
fid=fopen([opt.name '.mesh' '.epart.' num2str(Nsub)]);
A=textscan(fid,'%f');
fclose(fid);
ElPart=A{1,1};

fid=fopen([opt.name '.mesh' '.npart.' num2str(Nsub)]);
A=textscan(fid,'%f');
fclose(fid);
NDPart=A{1,1};

end

function opt = check_metis_options(opt)
    if ispc
        Default_metis_path='mpmetis.exe';
    else
        Default_metis_path='/mpmetis';
    end
    if isempty(opt)
        opt.metis_path=Default_metis_path;
        opt.name='tempMetis';
    end
    if isempty(opt.metis_path); 
        opt.metis_path=Default_metis_path; 
    end

end


function outraster = Equal_Area_Scott(raster, min_area, id, excludeid)
    % outraster = Equal_Area(inraster, min_area, Neq_areas, id)
    % This function divides a raster into equal area cells
    %
    % INPUTS:
    % raster      : raster to be devided. Each cell id in the raster is
    %                 treated as a different region
    % min_area      : The area of each equal area cell (Units are cells).
    % id            : Divides only the area associated with the particular id.
    %                 If empty then all areas are divided. 
    %
    % excludeid     : any area associated with ids in this list will not be divided
    %                   This is usefull to exclude no data cell.
    %
    % OUTPUTS
    % outraster     : The raster with the equal area cells.
    % 
    % Usage:
    % If the raster comes from ARCGIS layer then read it into the Matlab workspace using
    % [raster metadata]=readArcGisASCIIfile(asciiraster.asc);
    %
    % Next this function to create the equal area cells 
    %
    % and finaly convert them to ARCGIS raster ascii.asc file using 
    % createAsciiforRaster(filename,TAB,xl,yl,csz,nodata)
    % (The additional functions are part of the mSim software http://groundwater.ucdavis.edu/msim/)

    if isempty(id)
        id = unique(raster);
    end
    id = setdiff(id, excludeid);
    [nr nc] = size(raster);
    outraster = nan(nr,nc);

    for k = 1:length(id)
        inraster=raster==id(k);
        NumCells=sum(sum(inraster));
        NumCell_per_eq_area=round(min_area);
        Ncat=round(NumCells/NumCell_per_eq_area);
        tempraster=zeros(size(inraster));
        tempraster(~inraster)=nan;

        Strips=zeros(size(inraster));
        Nstrips=round(sqrt(Ncat));

        NumCell_strip=round(NumCells/Nstrips);
        cumsum_strips=cumsum(sum(inraster));
        istart=1;
        for i=1:Nstrips-1
            [~,iend]=min(abs(cumsum_strips-i*NumCell_strip));
            Strips(:,istart:iend)=i;
            istart=iend+1;
        end
        Strips(:,istart:end)=Nstrips;

        iCat=1;cntCells=0;
        up=true;
        for i=1:Nstrips
            fprintf('%g strip out of %g\n', [i Nstrips])
            %find the lower left unused pixel
            [Is Js]=ind2sub(size(Strips),find(Strips==i,1));
            if i~=Nstrips
                [Ie Je]=ind2sub(size(Strips),find(Strips==i+1,1));
            else
                Je=nc;
            end
            if up
                Is=1;Ie=nr;step=1;
            else
                Is=nr;Ie=1;step=-1;
            end
            for ii=Is:step:Ie
                for jj=Js:Je-1
                    if ~isnan(tempraster(ii,jj)) && tempraster(ii,jj)==0
                        tempraster(ii,jj)=iCat;
                        cntCells=cntCells+1;
                        if cntCells==NumCell_per_eq_area
                            iCat=iCat+1;
                            if iCat>Ncat
                                iCat=Ncat;
                            end
                            cntCells=0;
                        end
                    end
                end
            end
            up=~up;
        end
        outraster(inraster) = id(k)*1000+tempraster(inraster);
    end
end