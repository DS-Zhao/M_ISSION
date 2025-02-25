function [obs]=cut_slip(path_obs,Sites_Info)
%% Remove small arc segments and cycle jump data
% INPUT:
%     path_obs: R_rinex:Observation data storage path
%     Sites_Info: name and coordinate information of the stations
% OUTPUT:
%     obs: Modified observation data
% SAVE:
%      */OBS_cut/doy/sitedoy.mat:Processed OBS data
%% modified by Zhang P. et al., 2024/08
%% *****************************************
doys=Sites_Info.doy;
stations=Sites_Info.name;
current_doy=num2str(unique(doys));
path_cutobs = fullfile(fileparts(path_obs), 'raw_OBS_cut');
list_obs=dir([path_cutobs,'/',current_doy,'/*.mat']);
list_obs=now_files(stations,list_obs);
len=length(list_obs);
for i=1:len
    load([path_cutobs,'/',current_doy,'/',list_obs(i).name],'-mat');
    doy=list_obs(i).name(end-8:end-4);
    fields=fieldnames(obs);
    if ~isnan(find(strcmp(fields, 'GPSC1W' )))
        if ~(all(all(obs.GPSC1W==0)) || all(all(obs.GPSC2W==0)))
            obs=GPS_slip(obs);
        end
    end
    if ~isnan(find(strcmp(fields, 'GLOC1P' )))
        if ~(all(all(obs.GLOC1P==0)) || all(all(obs.GLOC2P==0)))
            obs=GLO_slip(obs);
        end
    end
    if ~isnan(find(strcmp(fields, 'BDSC2I')))
        if ~( all(all(obs.BDSC2I==0)) || all(all(obs.BDSC6I==0)))
            obs=BDS_slip(obs);
        end
    elseif ~isnan(find(strcmp(fields, 'BDSC1I')))
        if ~( all(all(obs.BDSC1I==0)) || all(all(obs.BDSC7I==0)))
            obs=BDS1_slip(obs);
        end
    end

    if ~isnan(find(strcmp(fields, 'GALC1C' )))
        if ~(all(all(obs.GALC1C==0)) || all(all(obs.GALC5Q==0)))
            obs=GAL_slip(obs);
        end
    elseif ~isnan(find(strcmp(fields, 'GALC1X' )))
        if ~(all(all(obs.GALC1X==0)) || all(all(obs.GALC5X==0)))
            obs=GALX_slip(obs);
        end
    end

    path_cutobs2 = fullfile(fileparts(path_obs), 'raw_OBS_cut',doy);
    filename1 = fullfile(path_cutobs2, list_obs(i).name);
    % Save data to file
    save(filename1, 'obs','-mat');
end
end
%%
function [obs]=GPS_slip(obs)
%%  GPS observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     obs: Modified observation data
%% --------------------------------------------------------------------------
%____delete incomplete epoch_____

size2=size(obs.GPSC1W,2);
GPS_f1=1575.42*10^6;                 %_________________________________unit:Hz
GPS_f2=1227.6*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(GPS_f1-GPS_f2);       %____________________wide lane wavelength
L6=lamda_w*(obs.GPSL1C-obs.GPSL2W)-(GPS_f1*obs.GPSC1W+GPS_f2*obs.GPSC2W)/(GPS_f1+GPS_f2); %__MW observable
Li=obs.GPSL1C-GPS_f1*obs.GPSL2W/GPS_f2;
Nw=-L6;                   %_____________________wide lane ambiguity
t=1.4;
for i=1:size2  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GPSC1W(k,i)=0;obs.GPSC2W(k,i)=0;obs.GPSL1C(k,i)=0;obs.GPSL2W(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GPSL1C(e,i)=0;obs.GPSL2W(e,i)=0;obs.GPSC1W(e,i)=0;obs.GPSC2W(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GPSC1W(k,i)=0;obs.GPSC2W(k,i)=0;obs.GPSL1C(k,i)=0;obs.GPSL2W(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<t*4*sigma(count)&&I1<t*0.28   %-------------------------no cycle slip
            % if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GPSC1W(k+1,i)=0;obs.GPSC2W(k+1,i)=0;
                        obs.GPSL1C(k+1,i)=0;obs.GPSL2W(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GPSC1W(l,i)=0;obs.GPSC2W(l,i)=0;
                            obs.GPSL1C(l,i)=0;obs.GPSL2W(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    %if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GPSC1W(l,i)=0;obs.GPSC2W(l,i)=0;
                            obs.GPSL1C(l,i)=0;obs.GPSL2W(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GPSC1W(k+1,i)=0;obs.GPSC2W(k+1,i)=0;
                        obs.GPSL1C(k+1,i)=0;obs.GPSL2W(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GPSC1W(l,i)=0;obs.GPSC2W(l,i)=0;
                            obs.GPSL1C(l,i)=0;obs.GPSL2W(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
end
end

%% ------------------------subfunction----------------------------
function [obs]=GLO_slip(obs)
%%  GLO observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     obs: Modified observation data
%% ------------------------------------------------------------------
%____delete incomplete epoch_____
size2=size(obs.GLOC1P,2);
for i=1:size2
    Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2,0,0];
    %      Fre=[1,-4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2];
    GLO_f1(i)=(1602+Fre(i)*0.5625)*10^6;                 %_________________________________unit:Hz
    GLO_f2(i)=(1246+Fre(i)*0.4375)*10^6;                  %_________________________________unit:Hz
    lamda_w(i)=299792458/(GLO_f1(i)-GLO_f2(i));       %____________________wide lane wavelength
    L6(:,i)=lamda_w(i)*(obs.GLOL1P(:,i)-obs.GLOL2P(:,i))-(GLO_f1(i)*obs.GLOC1P(:,i)+GLO_f2(i)*obs.GLOC2P(:,i))/(GLO_f1(i)+GLO_f2(i)); %__MW observable
    Li(:,i)=obs.GLOL1P(:,i)-GLO_f1(i)*obs.GLOL2P(:,i)/GLO_f2(i);
    Nw(:,i)=-L6(:,i);               %_____________________wide lane ambiguity
end
t=1.4;
for i=1:21  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GLOC1P(k,i)=0;obs.GLOC2P(k,i)=0;obs.GLOL1P(k,i)=0;obs.GLOL2P(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GLOL1P(e,i)=0;obs.GLOL2P(e,i)=0;obs.GLOC1P(e,i)=0;obs.GLOC2P(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GLOC1P(k,i)=0;obs.GLOC2P(k,i)=0;obs.GLOL1P(k,i)=0;obs.GLOL2P(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<t*4*sigma(count)&&I1<t*0.28
            % if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GLOC1P(k+1,i)=0;obs.GLOC2P(k+1,i)=0;
                        obs.GLOL1P(k+1,i)=0;obs.GLOL2P(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GLOC1P(l,i)=0;obs.GLOC2P(l,i)=0;
                            obs.GLOL1P(l,i)=0;obs.GLOL2P(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GLOC1P(l,i)=0;obs.GLOC2P(l,i)=0;
                            obs.GLOL1P(l,i)=0;obs.GLOL2P(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GLOC1P(k+1,i)=0;obs.GLOC2P(k+1,i)=0;
                        obs.GLOL1P(k+1,i)=0;obs.GLOL2P(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GLOC1P(l,i)=0;obs.GLOC2P(l,i)=0;
                            obs.GLOL1P(l,i)=0;obs.GLOL2P(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end

end
end

%% ---------------------------subfunction-------------------------
function [obs]=BDS_slip(obs)
%%  BDS observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     obs: Modified observation data
%% -----------------------------------------------------------------
%____delete incomplete epoch_____
size2=size(obs.BDSC2I,2);
BDS_f2=1561.098*10^6;                 %_________________________________unit:Hz
BDS_f6=1268.520*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(BDS_f2-BDS_f6);       %____________________wide lane wavelength
L6=lamda_w*(obs.BDSL2I-obs.BDSL6I)-(BDS_f2*obs.BDSC2I+BDS_f6*obs.BDSC6I)/(BDS_f2+BDS_f6); %__MW observable
Li=obs.BDSL2I-BDS_f2*obs.BDSL6I/BDS_f6;
Nw=-L6;                   %_____________________wide lane ambiguity
t=1.4;
for i=1:size2  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.BDSC2I(k,i)=0;obs.BDSC6I(k,i)=0;obs.BDSL2I(k,i)=0;obs.BDSL6I(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.BDSL2I(e,i)=0;obs.BDSL6I(e,i)=0;obs.BDSC2I(e,i)=0;obs.BDSC6I(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.BDSC2I(k,i)=0;obs.BDSC6I(k,i)=0;obs.BDSL2I(k,i)=0;obs.BDSL6I(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<t*4*sigma(count)&&I1<t*0.28
            % if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.BDSC2I(k+1,i)=0;obs.BDSC6I(k+1,i)=0;
                        obs.BDSL2I(k+1,i)=0;obs.BDSL6I(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.BDSC2I(l,i)=0;obs.BDSC6I(l,i)=0;
                            obs.BDSL2I(l,i)=0;obs.BDSL6I(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.BDSC2I(l,i)=0;obs.BDSC6I(l,i)=0;
                            obs.BDSL2I(l,i)=0;obs.BDSL6I(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.BDSC2I(k+1,i)=0;obs.BDSC6I(k+1,i)=0;
                        obs.BDSL2I(k+1,i)=0;obs.BDSL6I(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.BDSC2I(l,i)=0;obs.BDSC6I(l,i)=0;
                            obs.BDSL2I(l,i)=0;obs.BDSL6I(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
end
end
%%
function [obs]=BDS1_slip(obs)
%%  BDS observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     obs: Modified observation data
%% -----------------------------------------------------------------
size2=size(obs.BDSC1I,2);
BDS_f2=1561.098*10^6;               %_________________________________unit:Hz
BDS_f7=1207.140*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(BDS_f2-BDS_f7);       %____________________wide lane wavelength
L6=lamda_w*(obs.BDSL1I-obs.BDSL7I)-(BDS_f2*obs.BDSC1I+BDS_f7*obs.BDSC7I)/(BDS_f2+BDS_f7); %__MW observable
Li=obs.BDSL1I-BDS_f2*obs.BDSL7I/BDS_f7;
Nw=-L6;                   %_____________________wide lane ambiguity
t=1.4;
for i=1:size2  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.BDSC1I(k,i)=0;obs.BDSC7I(k,i)=0;obs.BDSL1I(k,i)=0;obs.BDSL7I(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.BDSL1I(e,i)=0;obs.BDSL7I(e,i)=0;obs.BDSC1I(e,i)=0;obs.BDSC7I(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.BDSC1I(k,i)=0;obs.BDSC7I(k,i)=0;obs.BDSL1I(k,i)=0;obs.BDSL7I(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<t*4*sigma(count)&&I1<t*0.28
            % if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.BDSC1I(k+1,i)=0;obs.BDSC7I(k+1,i)=0;
                        obs.BDSL1I(k+1,i)=0;obs.BDSL7I(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.BDSC1I(l,i)=0;obs.BDSC7I(l,i)=0;
                            obs.BDSL1I(l,i)=0;obs.BDSL7I(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.BDSC1I(l,i)=0;obs.BDSC7I(l,i)=0;
                            obs.BDSL1I(l,i)=0;obs.BDSL7I(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.BDSC1I(k+1,i)=0;obs.BDSC7I(k+1,i)=0;
                        obs.BDSL1I(k+1,i)=0;obs.BDSL7I(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.BDSC1I(l,i)=0;obs.BDSC7I(l,i)=0;
                            obs.BDSL1I(l,i)=0;obs.BDSL7I(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
end
end
%% --------------------------subfunction--------------------------
function [obs]=GAL_slip(obs)
%%  GAL observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     obs: Modified observation data
%% -----------------------------------------------------------------
%____delete incomplete epoch_____
size2=size(obs.GALC1C,2);
GAL_f1=1575.42*10^6;                 %_________________________________unit:Hz
GAL_f5=1176.45*10^6;                  %_________________________________unit:Hz
lamda_w=299792458/(GAL_f1-GAL_f5);       %____________________wide lane wavelength
L6=lamda_w*(obs.GALL1C-obs.GALL5Q)-(GAL_f1*obs.GALC1C+GAL_f5*obs.GALC5Q)/(GAL_f1+GAL_f5); %__MW observable
Li=obs.GALL1C-GAL_f1*obs.GALL5Q/GAL_f5;
Nw=-L6;                   %_____________________wide lane ambiguity
t=1.4;
for i=1:size2  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1C(k,i)=0;obs.GALC5Q(k,i)=0;obs.GALL1C(k,i)=0;obs.GALL5Q(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GALL1C(e,i)=0;obs.GALL5Q(e,i)=0;obs.GALC1C(e,i)=0;obs.GALC5Q(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1C(k,i)=0;obs.GALC5Q(k,i)=0;obs.GALL1C(k,i)=0;obs.GALL5Q(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<t*4*sigma(count)&&I1<t*0.28
            % if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1C(k+1,i)=0;obs.GALC5Q(k+1,i)=0;
                        obs.GALL1C(k+1,i)=0;obs.GALL5Q(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GALC1C(l,i)=0;obs.GALC5Q(l,i)=0;
                            obs.GALL1C(l,i)=0;obs.GALL5Q(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GALC1C(l,i)=0;obs.GALC5Q(l,i)=0;
                            obs.GALL1C(l,i)=0;obs.GALL5Q(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1C(k+1,i)=0;obs.GALC5Q(k+1,i)=0;
                        obs.GALL1C(k+1,i)=0;obs.GALL5Q(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GALC1C(l,i)=0;obs.GALC5Q(l,i)=0;
                            obs.GALL1C(l,i)=0;obs.GALL5Q(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
end
end

%% --------------------------subfunction----------------------------
function [obs]=GALX_slip(obs)
%%  GALX observations
% INPUT:
%      obs: struct of rinex files
% OUTPUT:
%     obs: Modified observation data
%% -----------------------------------------------------------------
%____delete incomplete epoch_____
size2=size(obs.GALC1X,2);
GALX_f1=1575.42*10^6;                 %_________________________________unit:Hz
GALX_f5=1176.45*10^6;                 %_________________________________unit:Hz
lamda_w=299792458/(GALX_f1-GALX_f5);       %____________________wide lane wavelength
L6=lamda_w*(obs.GALL1X-obs.GALL5X)-(GALX_f1*obs.GALC1X+GALX_f5*obs.GALC5X)/(GALX_f1+GALX_f5); %__MW observable
Li=obs.GALL1X-GALX_f1*obs.GALL5X/GALX_f5;
Nw=-L6;                   %_____________________wide lane ambiguity
t=1.4;
for i=1:size2  %i is PRN number
    %------divide arc---------------------------
    arc=Get_arc(L6(:,i));
    [arc_n,aaa]=size(arc);
    %----delete arc less than 10 epoches-------
    arc_d=[];
    for j=1:arc_n
        n_epoch=arc(j,2)-arc(j,1);
        if n_epoch<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1X(k,i)=0;obs.GALC5X(k,i)=0;obs.GALL1X(k,i)=0;obs.GALL5X(k,i)=0;
                L6(k,i)=0;Nw(k,i)=0;Li(k,i)=0;
            end
            arc_d=[arc_d,j];
        end
    end
    arc(arc_d,:)=[];
    %----mw detect cycle slip------------------
    [arc_n,aaa]=size(arc);
    j=1;
    while j<arc_n+1 %j is arc number
        %----first epoch check----------
        e=arc(j,1);
        while 1
            if e+1==arc(j,2) || e==arc(j,2)
                break;
            end
            fir=Nw(e,i);sec=Nw(e+1,i);thi=Nw(e+2,i);
            firl=Li(e,i);secl=Li(e+1,i);thil=Li(e+2,i);
            sub=abs(fir-sec);sub2=abs(sec-thi);
            subl=abs(firl-secl);subl2=abs(secl-thil);
            if sub>1||sub2>1||subl>1||subl2>1
                L6(e,i)=0;obs.GALL1X(e,i)=0;obs.GALL5X(e,i)=0;obs.GALC1X(e,i)=0;obs.GALC5X(e,i)=0;Nw(e,i)=0;Li(e,i)=0;
                e=e+1;
                arc(j,1)=e;
            else
                arc(j,1)=e;
                break;
            end
        end
        %----detect------------------
        if arc(j,2)-arc(j,1)<10
            for k=arc(j,1):arc(j,2)
                obs.GALC1X(k,i)=0;obs.GALC5X(k,i)=0;obs.GALL1X(k,i)=0;obs.GALL5X(k,i)=0;
                L6(k,i)=0;Nw(e,i)=0;Li(k,i)=0;
            end
            arc(j,:)=[];
            arc_n=arc_n-1;
            continue;
        end
        ave_N=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma2=linspace(0,0,arc(j,2)-arc(j,1)+1);
        sigma=linspace(0,0,arc(j,2)-arc(j,1)+1);
        ave_N(1)=Nw(arc(j,1),i);
        sigma2(1)=0;
        sigma(1)=0;
        count=2;
        for k=arc(j,1)+1:arc(j,2)-1 %----------------------check epoch k+1
            ave_N(count)=ave_N(count-1)+(Nw(k,i)-ave_N(count-1))/count;
            sigma2(count)=sigma2(count-1)+((Nw(k,i)-ave_N(count-1))^2-sigma2(count-1))/count;
            sigma(count)=sqrt(sigma2(count));
            T=abs(Nw(k+1,i)-ave_N(count));
            I1=abs(Li(k+1,i)-Li(k,i));
            if T<t*4*sigma(count)&&I1<t*0.28
            % if T<4*sigma(count)&&I1<0.28   %-------------------------no cycle slip
                count=count+1;
                continue;
            else
                if k+1==arc(j,2)            %---------------------arc end
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1X(k+1,i)=0;obs.GALC5X(k+1,i)=0;
                        obs.GALL1X(k+1,i)=0;obs.GALL5X(k+1,i)=0;Nw(k+1,i)=0;
                        Li(k,i)=0;arc(j,2)=k;
                    else                     %------delete scatter epoches
                        for l=arc(j,1):k+1
                            L6(l,i)=0;obs.GALC1X(l,i)=0;obs.GALC5X(l,i)=0;
                            obs.GALL1X(l,i)=0;obs.GALL5X(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,:)=[];
                        j=j-1;
                        arc_n=arc_n-1;
                    end
                    break;
                end
                I2=abs(Li(k+2,i)-Li(k+1,i));
                if abs(Nw(k+2,i)-Nw(k+1,i))<1&&I2<1%-----------cycle slip
                    if k+1-arc(j,1)>10
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+1,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else                    %------delete scatter epoches
                        for l=arc(j,1):k
                            L6(l,i)=0;obs.GALC1X(l,i)=0;obs.GALC5X(l,i)=0;
                            obs.GALL1X(l,i)=0;obs.GALL5X(l,i)=0;Nw(l,i)=0;
                            Li(k,i)=0;
                        end
                        arc(j,1)=k+1;
                        j=j-1;
                    end
                else                        %-----------------gross error
                    if k+1-arc(j,1)>10
                        L6(k+1,i)=0;obs.GALC1X(k+1,i)=0;obs.GALC5X(k+1,i)=0;
                        obs.GALL1X(k+1,i)=0;obs.GALL5X(k+1,i)=0;Nw(k+1,i)=0;Li(k,i)=0;
                        arc=[arc(1:j-1,:);[arc(j,1),k];[k+2,arc(j,2)];arc(j+1:arc_n,:)];
                        arc_n=arc_n+1;
                    else
                        for l=arc(j,1):k+1  %------delete scatter epoches
                            L6(l,i)=0;obs.GALC1X(l,i)=0;obs.GALC5X(l,i)=0;
                            obs.GALL1X(l,i)=0;obs.GALL5X(l,i)=0;Nw(l,i)=0;Li(k,i)=0;
                        end
                        arc(j,1)=k+2;
                        j=j-1;
                    end
                end
                break;
            end
        end
        j=j+1;
    end
end
end
%%
function arc = Get_arc(array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
len=length(array);
arc=[];
for i=1:len
    if i==len
        if array(i)~=0
            arc=[arc,i];
        end
        continue;
    end
    if i==1&&array(i)~=0
        arc=[arc,i];
    end
    if array(i)==0&&array(i+1)~=0
        arc=[arc,i+1];
        continue;
    end
    if array(i)~=0&&array(i+1)==0
        arc=[arc,i];
        continue;
    end
end
if len==0,return,end
arc=reshape(arc,2,[]);
arc=arc';
end
%%  Observation file corresponding file
function new_files=now_files(stations,files)
file_features = cell(length(files)-2, 1);  % Pre allocated cell array
for i = 1:length(files)
    file_name = files(i).name;
    file_features{i} = file_name(1:4);  % Extract the first four feature names
end
% Check if each feature name exists in stations
valid_files_indices = [];
for i = 1:length(file_features)
    if ismember(file_features{i}, stations)
        valid_files_indices = [valid_files_indices, i];
    end
end
new_files = files(valid_files_indices);
end
