function plot_Mulfiles(app,files_to_plot)
%% Draw multiple index files
%INPUT:
%     app: Software parameter input
%     Mfiles_to_plot:Storage paths for multiple index files
%                    files_to_plot{1}~{4}:ROTI,AATR,RMSAATR,SIGMAPHI
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
%Calculate the number of non-zero file paths
m = sum(cellfun(@(x) ~isequal(x, ''), files_to_plot));
%Initialize the current subgraph index
n = 1;
%Draw data for each non-zero path
if files_to_plot{1} ~= 0
    data = load(files_to_plot{1}, '-mat');
    name = cell2mat(fieldnames(data));
    namedata = data.(name);
    [validData,xValues]=time24h(namedata);
    [~, filename, ~] = fileparts(files_to_plot{1});
    ax1=subplot(m, 1, n,'parent',app.Panel_2);
    plot(ax1,xValues, validData);
    title(ax1,filename);
    xlim(ax1,[xValues(1) xValues(end)]);
    max_y = max(validData, [], 'all');
    min_y = min(validData, [], 'all');
    ylim(ax1,[min_y max_y]);
    if n==m
    xlabel(ax1,'Time (h)');
    end
    ylabel(ax1,'ROTI Value');
    n = n + 1;
end

if files_to_plot{2} ~= 0
    data = load(files_to_plot{2}, '-mat');
    name = cell2mat(fieldnames(data));
    namedata = data.(name);
    [validData,xValues]=time24h(namedata);
    [~, filename, ~] = fileparts(files_to_plot{2});
    ax2=subplot(m, 1, n,'parent',app.Panel_2);
    plot(ax2,xValues, validData);
    title(ax2,filename);
    xlim(ax2,[xValues(1) xValues(end)]);
    max_y = max(validData, [], 'all');
    min_y = min(validData, [], 'all');
    ylim(ax2,[min_y max_y]);
    if n==m
    xlabel(ax2,'Time (h)');
    end
    ylabel(ax2,'AATR Value');
    n = n + 1;
end


if files_to_plot{3} ~= 0
    data = load(files_to_plot{3}, '-mat');
    name = cell2mat(fieldnames(data));
    namedata = data.(name);
    [~, filename, ~] = fileparts(files_to_plot{3});
    ax3=subplot(m, 1, n,'parent',app.Panel_2);
    plot(ax3,namedata); 
    title(ax3,filename);
    xlim(ax3,[1 24]);
    max_y = max(namedata, [], 'all');
    min_y = min(namedata, [], 'all');
    ylim(ax3,[min_y max_y]);
    if n==m
    xlabel(ax3,'Time Point');
    end
    ylabel(ax3,'RMS AATR Value');
     n = n + 1;
end
if files_to_plot{4} ~= 0
    data = load(files_to_plot{4}, '-mat');
    name = cell2mat(fieldnames(data));
    namedata = data.(name);
    name2=fieldnames(namedata(1)); 
    if strcmp(name2(1),'L1')
       namedata=namedata.L1;
    else
       namedata=namedata.L2; 
    end
    [validData,xValues]=time24h(namedata);
    [~, filename, ~] = fileparts(files_to_plot{4});
    ax4=subplot(m, 1, n,'parent',app.Panel_2);
    plot(ax4,xValues, validData);
    title(ax4,filename);
    xlim(ax4,[xValues(1) xValues(end)]);
    max_y = max(validData, [], 'all');
    min_y = min(validData, [], 'all');
    ylim(ax4,[min_y max_y]);
    if n==m
    xlabel(ax4,'Time (h)');
    end
    ylabel(ax4,'sigmaphi Value');
end
end