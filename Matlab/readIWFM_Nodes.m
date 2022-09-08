function XY = readIWFM_Nodes(filename)

fid = fopen(filename, 'r');
read_header = true;
while 1
    try
        temp = fgetl(fid);
        if isempty(temp)
            continue
        end
        if strcmp(temp(1),'*')
            continue
        end
        if strcmp(temp(1),'C')
            continue
        end
        
        if read_header
            C = parseLine(temp);
            ND = str2double(C{1,1});
            C = parseLine(fgetl(fid));
            FACT = str2double(C{1,1});
            XY = nan(ND,3);
            read_header = false;
        else
            for ii = 1:ND
                if ii > 1
                    temp = fgetl(fid);
                end
                C = parseLine(temp);
                XY(ii,:) = [str2double(C{1,2}) str2double(C{1,3}) str2double(C{1,1})];
            end
            break;
        end
    catch
        break;
    end
    
end

fclose(fid);
end

function C = parseLine(line)
    C = strsplit(string(line),{'/', '\t', ' '});
    C = C(~cellfun('isempty',C));
end