function [xw, yw, well_found] = JitterWellLocation(xin,yin,pp,pos_acc, IRRXY,irr_thres, DOMXY, dom_thres, cnt_lim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Generate a location
xw = nan;
yw = nan;
well_found = false;
for k = 1:cnt_lim
    rxy = -pos_acc + (2*pos_acc)*rand(1,2);
    xtmp = xin + rxy(1);
    ytmp = yin + rxy(2);
    if ~isempty(pp)
        if ~pp.isinterior(xtmp, ytmp)
            continue;
        end
    end
    if ~isempty(IRRXY)
        dst = min(sqrt( (IRRXY(:,1) - xtmp).^2 + (IRRXY(:,1) - ytmp).^2 ));
        if dst < irr_thres
            continue
        end
    end

    if ~isempty(DOMXY)
        dst1 = min(sqrt( (DOMXY(:,1) - xtmp).^2 + (DOMXY(:,2) - ytmp).^2 ));
        if dst1 > dom_thres
            well_found = true;
            xw = xtmp;
            yw = ytmp;
            break;
        end
    else
        well_found = true;
        xw = xtmp;
        yw = ytmp;
        break;
    end

    
end

