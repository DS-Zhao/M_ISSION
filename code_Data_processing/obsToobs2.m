function [obss, coor] = obsToobs2(obs,inf)
%% PPPH reads version 2 observation files and converts them to the format of this software
% INPUT:
%     obs: Data obtained by PPPH
%     inf: PPPH section setting parameters
% OUTPUT:
%     obss: Converted observation values
%     coor:Receiver coordinates
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
    c = 299792458;
    
    % GPS
    if any(any(~isnan(obs.l1(:,1:32)))) && any(any(~isnan(obs.l2(:,1:1:32))))
        obss.GPSL1C = obs.l1(:,1:32) * (1575.42e6 / c);
        obss.GPSL2W = obs.l2(:,1:32) * (1227.6e6 / c);
        obss.GPSC1W = obs.p1(:,1:32);
        obss.GPSC2W = obs.p2(:,1:32);
    end
    
    % GLONASS
    if any(any(~isnan(obs.l1(:,33:56)))) && any(any(~isnan(obs.l2(:,33:56))))
        Fre = [1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
        GLO_f1 = (1602 + Fre' * 0.5625) * 1e6;
        GLO_f2 = (1246 + Fre' * 0.4375) * 1e6;
        obss.GLOL1P = obs.l1(:,33:56) .* (GLO_f1' / c);
        obss.GLOL2P = obs.l2(:,33:56) .* (GLO_f2' / c);
        obss.GLOC1P = obs.p1(:,33:56);
        obss.GLOC2P = obs.p2(:,33:56);
    end
    
    % Galileo
    if any(any(~isnan(obs.l1(:,59:88)))) && any(any(~isnan(obs.l2(:,59:88))))
                obss.GALL1X = obs.l1(:,59:88) * (1575.42e6 / c);
                obss.GALL5X = obs.l2(:,59:88) * (1227.6e6 / c);
                obss.GALC1X = obs.p1(:,59:88);
                obss.GALC5X = obs.p2(:,59:88); 
    end
    coor=inf.rec.pos;
    
end