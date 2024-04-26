function cft = NPSAT_lgnrmfit(urf,opt)
urf = urf(:);
if opt.usematlab
    
    ft = fittype( '(1/(x*b*sqrt(2*pi)))*exp((-(log(x)-a)^2)/(2*b^2))', 'independent', 'x', 'dependent', 'y' );
    optslocal = fitoptions( 'Method', 'NonlinearLeastSquares' );
    optslocal.Display = 'Off';
    optslocal.StartPoint = [0.499112701571756 0.336174051321482];

    X = 1:length(urf);
    [xData, yData] = prepareCurveData( X, urf );
    [fitresult, gof] = fit( xData, yData, ft, optslocal );
    cft = coeffvalues(fitresult);
else

    iter1 = 1;
    hasValidfid = false;
    while true
        fitname = [opt.fitname '_' num2str(opt.Eid) '_' num2str(opt.Sid) '_' num2str(round(10000000*rand)) '.dat'];
        fid1 = fopen(fitname, 'w');
        if fid1 >= 0
            hasValidfid = true;
            break;
        end
        iter1 = iter1 + 1;
        if iter1 > 20
            break;
        end
    end
    
    if hasValidfid
        fprintf(fid1,'%.10f\n', urf);
        fclose(fid1);
        [st, cm] = system([opt.exe_path ' ' fitname]);
        cft = textscan(strtrim(cm),'%f %f');
        delete(fitname);
    end
end

