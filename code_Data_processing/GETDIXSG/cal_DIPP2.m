function d=cal_DIPP2(F1,F2)
%% Calculate the distance between puncture points
% INPUT:
%     F1: Longitude and latitude of puncture point at station 1
%     F2: Longitude and latitude of puncture point at station 2
% OUTPUT:
%     d:Distance between puncture points
%% written by Zhang P. et al., 2024/08    
%% --------------------------------------------------------------

g=min(size(F1.B,2),size(F2.B,2));
d=NaN(2880,g); %The distance between two receivers corresponding to the IPP of the same satellite
R= 6828.137;
for i=1:g
    for j=1:2880-1
        if F1.L(j,i)~=0 && F2.L(j,i)~=0 && F1.B(j,i)~=0 && F2.B(j,i)~=0 && ~isnan(F1.B(j,i)) && ~isnan(F1.L(j,i)) && ~isnan(F2.L(j,i)) && ~isnan(F2.B(j,i))
            lon1 = (F1.L(j,i)) ; lat1 = (F1.B(j,i));  
            lon2 = (F2.L(j,i)) ; lat2 = (F2.B(j,i));  
            dlon = lat1 - lat2;
            dlat = lon1 - lon2;
            a = (sin(dlat./2)).^2 +  cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
            c = 2 * asin(sqrt(a));
            d(j,i)= R .* c;
        end
    end
end
end