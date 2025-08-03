function [data] = clk_jmp2(data)
%% Receiver Clock Jump Detection and Repair
% Function:
%     Detects and repairs millisecond-level jumps in the receiver clock. Such jumps,
%     often caused by receiver firmware, can severely impact precise positioning
%     and timing if not corrected.
%     The method is based on the ionosphere-free, geometry-free (IF-GF)
%     combination. By differentiating this combination between epochs, clock jumps
%     can be isolated and identified. A jump is confirmed if it is observed
%     simultaneously across multiple satellites with a magnitude close to an
%     integer number of milliseconds. The L1 and L2 carrier phase observations
%     are then repaired accordingly.
%
% INPUT:
%     data:       A structure containing GNSS data, which must include:
%                 - data.obs: A structure with dual-frequency observations
%                             (L1, L2, P1, P2).
%
% OUTPUT:
%     data:       The updated data structure with the following modifications:
%                 - data.obs.l1: L1 carrier phase observations are corrected
%                                following a detected clock jump.
%                 - data.obs.l2: L2 carrier phase observations are corrected
%                                following a detected clock jump.
%                 - data.jump:   A new field is added; a column vector that records
%                                the detected clock jump value for each epoch (in ms).
%
%% ---------------------------------------------------------------------
[arc]  = arc_dtr(data.obs);

sn = size(data.obs.st,2);

c   = 299792458;% m/s

pot  = NaN(size(data.obs.st,1),sn);
jump = zeros(size(data.obs.st,1),1);
lmt  = 1 - 50*10^-6; %millisecond w.r.t dif

for k=1:sn
    
    ark = arc{k};
    for t=1:size(ark,1)
        
        st = ark(t,1); fn = ark(t,2);
        if k<33 %for GPS
            ifp = i_free(data.obs.p1(st:fn,k),data.obs.p2(st:fn,k),0);
            ifl = i_free(data.obs.l1(st:fn,k),data.obs.l2(st:fn,k),0);
            df  = ifl - ifp;
        elseif k<59    %for GLONASS
            ifp = i_free(data.obs.p1(st:fn,k),data.obs.p2(st:fn,k),1);
            ifl = i_free(data.obs.l1(st:fn,k),data.obs.l2(st:fn,k),1);
            df  = ifl - ifp;
        elseif k<89
            ifp = i_free(data.obs.p1(st:fn,k),data.obs.p2(st:fn,k),2);
            ifl = i_free(data.obs.l1(st:fn,k),data.obs.l2(st:fn,k),2);
            df  = ifl - ifp;
        elseif k<106
            ifp = i_free(data.obs.p1(st:fn,k),data.obs.p2(st:fn,k),3);
            ifl = i_free(data.obs.l1(st:fn,k),data.obs.l2(st:fn,k),3);
            df  = ifl - ifp;   
        end
        
        dfi = diff(df)./(10^-3*c);
        if any(abs(dfi)>lmt)
            m = find(abs(dfi)>lmt);
            for i=1:size(m,1)
                if abs(round(dfi(m)) - dfi(m))<lmt
                    pot(st+m,k) = dfi(m);
                end
            end
        end
    end   
end


for i=1:size(pot,1)
    if any(~isnan(pot(i,:)))
        kern = find(~isnan(pot(i,:)));
        M = ((sum(pot(i,:),'omitnan')))/(length(kern));
        k2 = 10^-5; %ms
        if abs(M - round(M))<=k2
            jump(i,1) = round(M);
        end
    end
end

kern2 = find(jump~=0);
for k=1:sn
    ark = arc{k};
    for t=1:size(ark,1)
        st = ark(t,1); fn = ark(t,2);
        for ke = kern2'
            if (ke>st) && (ke<=fn)
                data.obs.l1(ke:fn,k) = data.obs.l1(ke:fn,k) + jump(ke,1)*(c*10^-3);
                data.obs.l2(ke:fn,k) = data.obs.l2(ke:fn,k) + jump(ke,1)*(c*10^-3);
            end
        end
    end
end
data.jump = jump;
end

