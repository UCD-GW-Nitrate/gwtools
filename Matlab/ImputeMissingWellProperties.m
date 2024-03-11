function well_values = ImputeMissingWellProperties(qv, dv, sv, pdf_QD, pdf_DS, pdf_QS, ...
    fitQD, fitDQ, fitDS, fitSD, fitQS, fitSQ, minmaxQ, minmaxD, minmaxS)
% ImputeMissingWellProperties: Fills in the missing properties of a well

bQ = isnan(qv);
bD = isnan(dv);
bS = isnan(sv);

minQ = minmaxQ(1);
maxQ = minmaxQ(2);
minD = minmaxD(1);
maxD = minmaxD(2);
minS = minmaxS(1);
maxS = minmaxS(2);

cnt = 1;

cnt_lim = 2000;
lim_reach = true;
if ~bQ && ~bD && ~bS
    well_values = [qv dv sv];
elseif bQ && bD && bS % All three are missing
    while cnt < cnt_lim
        qr = minQ + (maxQ - minQ)*rand;
        dr = minD + (maxD - minD)*rand;
        sr = minS + (maxS - minS)*rand;
        r1 = rand;
        r2 = rand;
        r3 = rand;
        if r1 < pdf_QS.F(qr,sr) && r2 < pdf_DS.F(dr,sr) && r3 < pdf_QD.F(qr,dr)
            well_values = [10^qr 10^dr 10^sr];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        dr = minD + (maxD - minD)*rand;
        well_values = [10^fitDQ(dr) 10^dr 10^fitDS(dr)];
    end
elseif bQ && bD && ~bS % Two missing
    while cnt < cnt_lim
        qr = minQ + (maxQ - minQ)*rand;
        dr = minD + (maxD - minD)*rand;
        r1 = rand;
        r2 = rand;
        r3 = rand;
        if r1 < pdf_QS.F(qr, log10(sv)) && r2 < pdf_DS.F(dr,log10(sv)) && r3 < pdf_QD.F(qr,dr)
            well_values = [10^qr 10^dr sv];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        well_values = [10^fitSQ(log10(sv)) 10^fitSD(log10(sv)) sv];
    end
elseif bQ && ~bD && bS % Two missing
    while cnt < cnt_lim
        qr = minQ + (maxQ - minQ)*rand;
        sr = minS + (maxS - minS)*rand;
        r1 = rand;
        r2 = rand;
        r3 = rand;
        if r1 < pdf_QS.F(qr,sr) && r2 < pdf_DS.F(log10(dv),sr) && r3 < pdf_QD.F(qr,log10(dv))
            well_values = [10^qr dv 10^sr];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        well_values = [10^fitDQ(log10(dv)) dv 10^fitDS(log10(dv))];
    end
elseif ~bQ && bD && bS % Two missing
    while cnt < cnt_lim
        dr = minD + (maxD - minD)*rand;
        sr = minS + (maxS - minS)*rand;
        r1 = rand;
        r2 = rand;
        r3 = rand;
        if r1 < pdf_QS.F(log10(qv),sr) && r2 < pdf_DS.F(dr,sr) && r3 < pdf_QD.F(log10(qv), dr)
            well_values = [qv 10^dr 10^sr];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        well_values = [qv 10^fitQD(log10(qv)) 10^fitQS(log10(qv))];
    end
elseif bQ && ~bD && ~bS % One missing
    while cnt < cnt_lim
        qr = minQ + (maxQ - minQ)*rand;
        r1 = rand;
        r2 = rand;
        if r1 < pdf_QS.F(qr,log10(sv)) && r2 < pdf_QD.F(qr,log10(dv))
            well_values = [10^qr dv sv];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        well_values = [10^fitDQ(log10(dv)) dv sv];
    end
elseif ~bQ && bD && ~bS % One missing
    while cnt < cnt_lim
        dr = minD + (maxD - minD)*rand;
        r1 = rand;
        r2 = rand;
        if r1 < pdf_DS.F(dr,log10(sv)) && r2 < pdf_QD.F(log10(qv),dr)
            well_values = [qv 10^dr sv];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        well_values = [qv 10^fitSD(log10(sv)) sv];
    end
elseif ~bQ && ~bD && bS % One missing
    while cnt < cnt_lim
        sr = minS + (maxS - minS)*rand;
        r1 = rand;
        r2 = rand;
        if r1 < pdf_QS.F(log10(qv),sr) && r2 < pdf_DS.F(log10(dv),sr)
            well_values = [qv dv 10^sr];
            lim_reach = false;
            break
        else
            cnt = cnt + 1;
        end
    end
    if lim_reach
        well_values = [qv dv 10^fitDS(log10(dv))];
    end
else
    warning('There should not be any other option')
    well_values = nan(1,3);
end