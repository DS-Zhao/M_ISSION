function [dif,B,L]=Get_dif(sx,sy,sz,E,A)
%% Calculate the latitude and longitude of the measuring point
% INPUT:
%      sx,sy,sz:Receiver three coordinates
%      E:satellite elevation
%      A:satellite azimuth
% OUTPUT:
%     B: latitude of Puncture point
%     L: Longitude of Puncture point 
%     dif:Distance between puncture points
%% written by Zhang P. et al., 2024/08    
%% -------------------------------------------------------------- 
        [sb,sl]=XYZtoBLH(sx,sy,sz);
        [m,n]=size(E);
        IPPz=NaN(m,n);t=NaN(m,n);b=NaN(m,n);ss=NaN(m,n);B=NaN(m,n);L=NaN(m,n);
        for i=1:n
            for j=1:m
                if E(j,i)==0
                    continue
                else
                    IPPz(j,i)=asin(6371000*sin(pi/2-E(j,i))/(6371000+450000));
                    t(j,i)=pi/2-E(j,i)-IPPz(j,i);
                    b(j,i)=asin(sin(sb).*cos(t(j,i))+cos(sb).*sin(t(j,i)).*cos(A(j,i)));
                    ss(j,i)=sl+asin(sin(t(j,i)).*sin(A(j,i))./cos(b(j,i)));
                    B(j,i)=b(j,i);
                    L(j,i)=ss(j,i);
                end
            end
        end
        dif=cal_IPP2(B,L);
end