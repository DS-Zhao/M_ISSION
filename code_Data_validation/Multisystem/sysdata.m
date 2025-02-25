function GREC=sysdata(Mfiles_to_plot)
%% Identify which system data exists in multi system data
%INPUT:
%     Mfiles_to_plot: Multi system data file path
%OUTPUT:
%     GREC:Identify existing system data
%% written by Zhang P. et al., 2024/11
%% -----------------------------------------------------------------------
if isfile(Mfiles_to_plot)
    GREC=[0 0 0 0];
    data = load(Mfiles_to_plot, '-mat');
    name = cell2mat(fieldnames(data));
    namedata = data.(name);
    fields = fieldnames(namedata);
    m = size(fields, 1);
    for i = 1:m
        GRECname=fields{i};
        if strcmp(GRECname(1:3),'GPS')
            GREC(1)=1;
        elseif strcmp(GRECname(1:3),'GLO')
            GREC(2)=1;
        elseif strcmp(GRECname(1:3),'GAL')
            GREC(3)=1;
        elseif strcmp(GRECname(1:3),'BDS')
            GREC(4)=1;
        end
    end
end