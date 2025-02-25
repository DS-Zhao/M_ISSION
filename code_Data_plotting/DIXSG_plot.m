function DIXSG_plot(app,DIXSGpath)
%% Draw DIXSG data file
%Draw hourly regional maps or histograms
%INPUT:
%     app: Software parameter input
%     DIXSGpath: The path of the DIXSGpath data file
%% written by Zhang P. et al., 2024/08
%% -----------------------------------------------------------------------
if exist(DIXSGpath,"file")
    load(DIXSGpath,'-mat');
            choice = questdlg('Select Shape:', ...
                'Select Shape', ...
                'A certain hour', 'DIXSG', 'cancel', 'cancel');
            
            %Handling user choices
            switch choice
                case 'A certain hour'
                    code=1;
                    prompt = {'enter the number of hours'};
                    dlgTitle = '';
                    numLines = 1;
                    defaultInput = {'1'};
                    userInput = inputdlg(prompt, dlgTitle, numLines, defaultInput);
                    if isempty(userInput)
                        return;
                    end 
                    num= str2double(userInput{1});
                    ax=figure;
                    %Draw a regional map
                     hourDIXSG(MBL(5),LL,num,MBL(2),MBL(4),MBL(1),MBL(3),ax);
                case 'DIXSG'
                    %Draw a histogram
                    DIXSG_P(app,aDIXSG); 
                case 'cancel'
                    s=sprintf('%s','cancel');
                    app.EditField.Value=s;
            end
end
end