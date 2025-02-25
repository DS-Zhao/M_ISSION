function [plotdata,xValues]=timedata(plottime,data)
%% Extract data corresponding to time and determine time labels
%INPUT:
%     plottim: The time interval drawn
%     data:Data of the corresponding system
%OUTPUT:
%     plotdata:Selected data to be drawn
%     xValues:The x-axis time label corresponding to the valid data
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
n=size(data,1);
h1=plottime(1);m1=plottime(2);h2=plottime(3);m2=plottime(4);
if h1<=h2 && h1*60+m1<h2*60+m2
    if h1==0 && h2==23 && m1==0 && m2==60
        plotdata=data;
        [plotdata,xValues]=time24h(plotdata);
    else
        datahead=n/24*(h1+m1/60)+1;
        dataend=n/24*(h2+m2/60);
        plotdata=data(datahead:dataend,:);
        xValues = (datahead:dataend)' * 24 / n;
    end
else
    plotdata=data;
    [plotdata,xValues]=time24h(plotdata);
end
end