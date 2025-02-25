function ROTI_plot(app, ROTIpath, file,allnum)
%% Draw a single ROTI file
%INPUT:
%     app: Software parameter input
%     ROTIpath: The path of the generated single ROTI data file
%     file：file name
%     allnum：Satellite number to be drawn
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
loadroti = load(ROTIpath, '-mat');
name = cell2mat(fieldnames(loadroti));
namedata=loadroti.(name);
if ~strcmp(allnum, 'all')
    allnum = strrep(allnum, '[', ''); % Remove the left bracket
    allnum = strrep(allnum, ']', ''); % Remove the right bracket
    allnum = strsplit(allnum, ','); % Separate strings based on commas
    a = [];
    for i = 1:length(allnum)
        if contains(allnum{i}, ' ')
            range = str2num(allnum{i}); % Convert a string containing spaces into an array
            a = [a, range(1):range(2)]; % Add numbers within the range to the array
        else
            a = [a, str2num(allnum{i})]; % Add a single number to an array
        end
    end
    a = unique(a); % Remove duplicate values
    a=sort(a);
    namedata=namedata(:,a);
end
[namedata,xValues]=time24h(namedata);



if size(namedata,2)<11
    % Draw each column of data
    hold(app.UIAxes, 'on'); 
    % Draw a chart in UIAxes
    colors = parula(size(namedata, 2)); 
    for col = 1:size(namedata, 2)
        plot(app.UIAxes, xValues, namedata(:, col), 'Color', colors(col,:), 'DisplayName', ['sat', num2str(a(col))]);
        legendEntries{col} = ['sat ', num2str(a(col))]; % Save column numbers for legend
    end
    legend(app.UIAxes, legendEntries, 'Location', 'Best');
    hold(app.UIAxes, 'off');
else
    plot(app.UIAxes, xValues, namedata);  
end
    %Set chart title and axis labels
    title(app.UIAxes, file);
    xlim(app.UIAxes, [xValues(1) xValues(end)]);
    max_y = max(namedata, [], 'all'); 
    min_y = min(namedata, [], 'all');
    ylim(app.UIAxes,[min_y-0.1*abs(min_y) max_y+0.1*abs(max_y)]);
    xlabel(app.UIAxes, 'Time(h)');
    ylabel(app.UIAxes, 'ROTI Value');
    xticks(app.UIAxes, 'auto');
    xticklabels(app.UIAxes, 'auto');
    %Add interactive zoom in and pan functions
    zoom(app.UIAxes, 'on');
end
