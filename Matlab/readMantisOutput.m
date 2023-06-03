function [btc, tf] = readMantisOutput(filename)
% readMantisOutput Reads the outfile of the TestClient program.
% Returns an Nbtc x Nyrs matrix where Nbtc is the number of well
% breakthrough curves and Nyrs is the number of simulation years.
fid = fopen(filename,'r');
temp = textscan(fid, '%f %f',1);
Nbtc = temp{1};
Nyrs = temp{2};

if Nbtc == 0
    tf = false;
    temp = fgetl(fid);
    display(temp);
    btc = [];
else
    tf = true;
    temp = textscan(fid, '%f',Nbtc*Nyrs);
    btc = reshape(temp{1},Nyrs,Nbtc)';
end
fclose(fid);
end

