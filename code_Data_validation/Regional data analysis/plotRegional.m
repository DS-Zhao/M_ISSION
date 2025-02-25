function plotRegional(Regionalpath,hour)
%% Draw DIXSG and ROTI data for the region
%INPUT:
%     Regionalpath: DIXSG and ROTI data paths in the region
%     hour:Draw the data for the hour
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
n = sum(~cellfun(@isempty, Regionalpath));
plotdixsg_path=Regionalpath{1};
plotRegROTI_path=Regionalpath{2};
i=1;

if exist(plotdixsg_path,"file")
    ax1=subplot(n,1,i);
    load(plotdixsg_path,'-mat');
    hourDIXSG(MBL(5),LL,hour,MBL(2),MBL(4),MBL(1),MBL(3),ax1);
    i=i+1;
end
if exist(plotRegROTI_path,"file")
    ax2=subplot(n,1,i);
    load(plotRegROTI_path,'-mat');
    hourROTI(ROTIgrid, hour, lonlat(3), lonlat(1), lonlat(4), lonlat(2),ax2);
end
end