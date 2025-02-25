function SIGMAPHI_plot(app,full_path,file,allnum)
%% Draw a single SIGMAPHI file
%Please refer to the ROTI_plot function for detailed annotations
%INPUT:
%     app: Software parameter input
%     full_path: The path of the generated single SIGMAPHI data file
%     file：file name
%     allnum：Satellite number to be drawn
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
sigmaphis=load(full_path);
sigmaphis = struct2cell(sigmaphis);
sigmaphis = sigmaphis{1};
if ismember('L1', fieldnames(sigmaphis))
    namedata=sigmaphis.L1;
else
    namedata=sigmaphis.L2;
end
if ~strcmp(allnum, 'all')
allnum = strrep(allnum, '[', ''); 
allnum = strrep(allnum, ']', ''); 
allnum = strsplit(allnum, ','); 
a = [];
for i = 1:length(allnum)
    if contains(allnum{i}, ' ')
        range = str2num(allnum{i}); 
        a = [a, range(1):range(2)]; 
    else
        a = [a, str2num(allnum{i})]; 
    end
end
a = unique(a); 
a=sort(a);
namedata=namedata(:,a);
end
[namedata,xValues]=time24h(namedata); 
if size(namedata,2)<11
    
    hold(app.UIAxes, 'on'); 
   
    colors = parula(size(namedata, 2)); 
    for col = 1:size(namedata, 2)
        plot(app.UIAxes, xValues, namedata(:, col), 'Color', colors(col,:), 'DisplayName', ['sat', num2str(a(col))]);
        legendEntries{col} = ['sat ', num2str(a(col))]; 
    end
    legend(app.UIAxes, legendEntries, 'Location', 'Best');
    hold(app.UIAxes, 'off'); 
else
plot(app.UIAxes, xValues, namedata);  
end

title(app.UIAxes, file);
xlim(app.UIAxes, [xValues(1) xValues(end)]);
max_y = max(namedata, [], 'all'); 
min_y = min(namedata, [], 'all');
ylim(app.UIAxes,[min_y-0.1*abs(min_y) max_y+0.1*abs(max_y)]);
xlabel(app.UIAxes, 'Time (h)');
ylabel(app.UIAxes, 'SIGMAPHI Value');
xticks(app.UIAxes, 'auto');
xticklabels(app.UIAxes, 'auto');

zoom(app.UIAxes, 'on');
end