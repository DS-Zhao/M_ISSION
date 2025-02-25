function DIXSG_P(app, DIXSG)
%% Draw a 24-hour histogram
%INPUT:
%     app: app: Software parameter input
%     DIXSG:24-hour DIXSG value
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
% Draw a chart in UIAxes
bar(app.UIAxes, 0:length(DIXSG)-1, DIXSG);
%Set chart title and axis labels
title(app.UIAxes, 'DIXSG');
xlabel(app.UIAxes, 'Time (h)');
ylabel(app.UIAxes, 'DIXSG');
xlim(app.UIAxes, [-1 length(DIXSG)+0.5]);
max_y = max(DIXSG); %Get maximum value
ylim(app.UIAxes,[0 max_y*1.1]);
xticks(app.UIAxes, 0:24);
xticklabels(app.UIAxes, strcat(string(1:24), 'h'));
zoom(app.UIAxes, 'on');
end