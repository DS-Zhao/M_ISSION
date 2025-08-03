function [data] = data_hand(files,options)
%% Read and integrate all input files into a single data structure
% Function:
%     This function serves as the main data loading entry point. It reads 
%     various GNSS data files such as Observation (OBS), Precise Ephemeris (SP3), 
%     Precise Clock (CLK), and Antenna Exchange (ATX) files based on the 
%     provided paths and configuration options. All read data is then 
%     consolidated into a single output structure 'data' for subsequent processing.
%
% INPUT:
%     files:    A structure containing the file paths for all necessary input data.
%               e.g., files.obs, files.orbit, files.clock, files.anten
%     options:  A structure containing various configuration options for data processing.
%               e.g., options.clock specifies the source of clock corrections.
%
% OUTPUT:
%     data:     A structure containing all the processed data, with the following main fields:
%               - inf:    Basic information extracted from file headers.
%               - obs:    Observation data.
%               - sat:    Satellite orbit data.
%               - clk:    Satellite clock correction data.
%               - atx:    Antenna Phase Center (APC) correction data.
%               - opt:    The input configuration options structure.
%               - files:  The input file paths structure.
%% ---------------------------------------------------------------------
narginchk(2,2)

[obs,inf] = read_obsf(files,options);

[sat,inf] = r_sp3(files.orbit,options,inf);

options.entrp = 0; 
if ~isempty(files.orbitb) && ~isempty(files.orbita)
    
    [satb,~] = r_sp3(files.orbitb,options,inf);
    
    [sata,~] = r_sp3(files.orbita,options,inf);

    if (size(satb.sp3,1)==size(sata.sp3,1))&&(size(satb.sp3,1)==size(sat.sp3,1))...
            &&(size(sat.sp3,1)==size(sata.sp3,1))
        sat.sp3 = vertcat(satb.sp3(end-4:end,:,:),...
                          sat.sp3(:,:,:),...
                          sata.sp3(1:5,:,:));
    end
    options.entrp = 1;
end

if strcmp(options.clock,'Clock File')
    [clk] = r_clck(files.clock,options);
    inf.time.clkint = options.clck_int;
elseif strcmp(options.clock,'Sp3 File')
    [clk] = sat.sp3(:,4,:);
    inf.time.clkint = inf.time.sp3int;
end

[atx] = r_antx(files.anten,inf,options);

data.inf  = inf;
data.robs = obs;
data.obs  = obs;
data.sat  = sat;
data.clk  = clk;
data.atx  = atx;
data.opt  = options;
data.files = files;

% save('data.mat', 'data', '-mat');
end