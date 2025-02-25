function dif=cal_IPP2(BB,LL)
%% Calculate the distance between puncture points
% haversine equation
% INPUT:
%     BB: latitude of Puncture point
%     LL: Longitude of puncture point
% OUTPUT:
%     dif:Distance between puncture points
%% written by Zhang X. et al., 2024/08    
%% --------------------------------------------------------------
g=size(BB,2);
dif=NaN(2880,g);
R=6828.1370;
for i=1:g
    for j=1:2880-1
        if ~isnan(LL(j,i)) && ~isnan(LL(j+1,i)) && ~isnan(BB(j,i)) && ~isnan(BB(j+1,i)) && LL(j,i) ~= 0 && LL(j+1,i) ~= 0 && BB(j,i) ~= 0 && BB(j+1,i) ~= 0
            lon1 = (LL(j+1,i)) ; lat1 = (BB(j+1,i));  %p1
            lon2 = (LL(j,i))   ; lat2 = (BB(j,i));   %p2
            dlon = lat1 - lat2;
            dlat = lon1 - lon2;
            a = (sin(dlat./2)).^2 +  cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
            c = 2 * asin(sqrt(a));
            dif(j,i)= R .* c;
        end
    end
end