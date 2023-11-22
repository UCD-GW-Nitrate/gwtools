function urfs = expandNPSATURF(m,s,Exit,NsimYears)
% urf = expandNPSATURF(m,s,Exit,NsimYears)
%
% Calculate the urf as a function of mean m and standard deviation s. 
% If m is -1 or -2 it uses a standard form for a URF of 1 and 2
% years respectively. 
% It calculates the URF only for those that the Exit is 1

OneYear = [0.48443252, 0.41642340, 0.08307405, 0.01338364, 0.00219913, 0.00039049, 0.00007584, 0.00001608, 0.00000370, 0.00000092, 0.00000023];
TwoYears = [0.00095171, 0.19513243, 0.41957050, 0.25244126, 0.09333424, 0.02803643, 0.00771540, 0.00206063, 0.00055016, 0.00014916, 0.00004141, 0.00001182, 0.00000347, 0.00000106, 0.00000032];

urfs = zeros(length(m),NsimYears);
for ii = 1:length(m)
    if Exit(ii) == 1
        if m(ii) == -1
            urfs(ii,:) = [OneYear zeros(1,NsimYears-length(OneYear))];
        elseif m(ii) == -2
            urfs(ii,:) = [TwoYears zeros(1,NsimYears-length(TwoYears))];
        else
            urfs(ii,:) = lognpdf(1:NsimYears,m(ii), s(ii));
        end
    end
end
end