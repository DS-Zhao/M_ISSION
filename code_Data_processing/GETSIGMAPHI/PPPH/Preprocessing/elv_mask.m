function [data] = elv_mask(data,options)
%% Apply Elevation Mask and Calculate Azimuth/Elevation Angles
% Function:
%     Calculates the azimuth and elevation angles for each satellite-to-receiver
%     line-of-sight. It then applies an elevation mask based on a user-defined
%     minimum angle, marking observations from satellites below this threshold
%     as invalid by setting their status to 0 in `data.obs.st`.
%
% INPUT:
%     data:       A structure containing GNSS data, which must include:
%                 - data.obs.st: The observation status matrix.
%                 - data.psat:   Satellite positions in the ECEF frame.
%                 - data.inf.rec.pos: Receiver position in the ECEF frame.
%     options:    A configuration structure that must contain the field:
%                 - options.elvangle: The minimum elevation angle threshold (in degrees).
%
% OUTPUT:
%     data:       The updated data structure with the following modifications:
%                 - data.obs.st: The status matrix is updated, flagging low-elevation
%                                satellites with a status of 0.
%                 - data.obs.elv: A new field storing the elevation angles for all
%                                 satellites (in degrees).
%                 - data.obs.azm: A new field storing the azimuth angles for all
%                                 satellites (in degrees).
%% ---------------------------------------------------------------------
en = size(data.obs.st,1);
sn = size(data.obs.st,2);

r_xyz = data.inf.rec.pos;
elv = NaN(en,sn);
azm = NaN(en,sn);
for i=1:en
    for k=1:sn
        if data.obs.st(i,k) == 1
            
            s_xyz = data.psat(i,1:3,k);
            
            [az,elev] = local(r_xyz,s_xyz,1);
            elv(i,k) = elev;
            azm(i,k) = az;
            if elev<options.elvangle
                data.obs.st(i,k) = 0;
            end
        end
    end
end

data.obs.elv = elv;
data.obs.azm = azm;
end

