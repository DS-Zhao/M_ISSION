function plot_Mfiles(app,Mfiles_to_plot,plottime,plotsys)
%% Draw multi system data: GPS,GAL,GLO,BDS
% Should choose files such as M.ROTI.mat, M.AATR, and M.RMSAATR
%INPUT:
%     app: Software parameter input
%     Mfiles_to_plot:Multi system file path
%     plottime:The time interval drawn
%     plotsys:The system that needs to be drawn
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
if isfile(Mfiles_to_plot)
    data = load(Mfiles_to_plot, '-mat');
    name = cell2mat(fieldnames(data));
    namedata = data.(name);
    fields = fieldnames(namedata);
    m = size(fields, 1);
    for i = 1:m
        GRECname = fields{i};
        if strcmp(GRECname(1:3),'GPS') && plotsys(1)==1    
            plotsysdata.(GRECname)=namedata.(GRECname);
        elseif strcmp(GRECname(1:3),'GLO') && plotsys(2)==1
            plotsysdata.(GRECname)=namedata.(GRECname);
        elseif strcmp(GRECname(1:3),'GAL') && plotsys(3)==1
            plotsysdata.(GRECname)=namedata.(GRECname);
        elseif strcmp(GRECname(1:3),'BDS') && plotsys(4)==1
            plotsysdata.(GRECname)=namedata.(GRECname);
        end
    end

    fields = fieldnames(plotsysdata);
    n = size(fields, 1);
    xValues=1:24;
    for i = 1:n
        %Extract data and draw subgraphs
        ax=subplot(n, 1, i,'parent',app.Panel);
        structname = fields{i};
        plotdata = plotsysdata.(structname);
        if ~strcmp(name,'M_RAATR')
        [plotdata,xValues]=timedata(plottime,plotdata);
        end
        plot(ax,xValues, plotdata);
        title(ax,structname);

        %Adjust the x and y axis range of the subgraph
        xlim(ax,[xValues(1) xValues(end)]);
        max_y = max(plotdata, [], 'all');
        min_y = min(plotdata, [], 'all');
        ylim(ax,[min_y, max_y]);

       if i==n
        xlabel(ax,'Time (h)');
       end
        ylabel(ax,'Value');
    end
end
end