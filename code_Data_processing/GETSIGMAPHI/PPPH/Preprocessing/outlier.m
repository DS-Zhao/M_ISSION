function [data] = outlier(data)
%% Outlier Detection and Removal
% Function:
%     This function detects and removes outliers (gross errors) from dual-frequency
%     GNSS observation data. The method is based on the Melbourne-Wübbena (MW)
%     combination, which should remain stable or vary smoothly over a continuous
%     arc in the absence of cycle slips and outliers.
%
%     Algorithm:
%     1. Process each satellite's data on an arc-by-arc basis.
%     2. For each continuous observation arc, fit a quadratic polynomial to the
%        time series of the MW combination values.
%     3. Calculate the residuals (difference between observed and fitted values)
%        for each epoch.
%     4. Iteratively identify and flag the epoch with the largest residual as an
%        outlier. This process repeats until all residuals are below a threshold
%        or the overall goodness-of-fit is acceptable.
%     5. The status flag `data.obs.st` for the detected outlier is set to 0.
%
% INPUT:
%     data:       A structure containing GNSS data, which must include:
%                 - data.obs: A structure with L1,L2,P1,P2 observations and a
%                             status flag `st`.
%
% OUTPUT:
%     data:       The updated data structure. This function modifies the
%                 `data.obs.st` field in place, setting the status flag to 0
%                 for any epoch identified as an outlier, effectively removing
%                 it from subsequent processing.
%
%% ---------------------------------------------------------------------
c  = 299792458; % m/s
[freq,~] = frequencies;

[arc]  = arc_dtr(data.obs);


sn = size(data.obs.st,2);


for k=1:sn
    
    f1 = freq(k,1); f2 = freq(k,2);
    
    lwl = (data.obs.l1(:,k).*f1 - data.obs.l2(:,k).*f2)./(f1-f2);
    pnl = (data.obs.p1(:,k).*f1 + data.obs.p2(:,k).*f2)./(f1+f2);
    lamwl = c/(f1 - f2);
    nwl = (lwl - pnl)./(lamwl);
    
    ark = arc{k};
    
    for n=1:size(ark,1)
        st = ark(n,1);
        fn = ark(n,2);
        while 1
            t = find(data.obs.st(:,k)==1);
            t = t(t(:)>=st & t(:)<=fn);
            L = nwl(t);
            ran = size(t,1);
            if(ran<3)
                data.obs.st(st:fn,k) = 0;
                break;
            end
            A = [(t).^2 t ones(ran,1)];
            
            X = A\L;
            V = L - A*X;
            
            rmse = sqrt(sum(V.^2)/ran-3);
            det = find(abs(V)>(4*0.6));
            if rmse>0.6 && any(det)
                
                del = t(abs(V)==max(abs(V)));
                data.obs.st(del,k) = 0;
            else
                break
            end
        end
    end
end
end

