function [final_xyz, success_flag, message] = glab_ppp_solver(data)
%GLAB_PPP_SOLVER Solves for position with a prioritized coordinate source chain.
%
%   Determines coordinates with the following priority:
%   1. Looks up the station in a custom coordinates text file.
%   2. If not found, runs gLAB PPP (with robust parsing and antenna error fallback).
%   3. If gLAB fails, uses the provided approximate coordinates.
%
%   SYNTAX:
%   [final_xyz, success_flag, message] = glab_ppp_solver(data)
%
%   INPUTS:
%   data: A structure containing:
%         .files.rinex: Path to RINEX observation file.
%         .files.orbit: Path to SP3 orbit file.
%         .files.clock: Path to CLK clock file.
%         .files.anten: Path to ANTEX antenna file.
%         .inf.rec.pos: (1x3 double) Approximate ECEF [X,Y,Z] as fallback.
%
%   OUTPUTS:
%   final_xyz:    (1x3 double) The final determined ECEF position.
%   success_flag: (logical) `true` on success, `false` on failure.
%   message:      (string) A detailed status message.
%
% Author: AI Assistant (based on user's evolving requirements)
% Date: 2025-06-30 (Robust Parsing Logic Integration)

% --- Setup paths and inputs from data structure ---
input_files.rinex_obs = data.files.rinex;
input_files.sp3       = data.files.orbit;
input_files.clk       = data.files.clock;
input_files.ant       = data.files.anten;
approx_xyz = data.inf.rec.pos;
custom_file=data.files.XYZ; 
removeRinexSplice(input_files.rinex_obs);
glab_paths.executable = 'gLAB.exe';
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
glab_paths.working_dir = script_dir;

if ~isnumeric(approx_xyz) || numel(approx_xyz) ~= 3
    final_xyz = [NaN, NaN, NaN];
    success_flag = false;
    return;
end

final_xyz = approx_xyz;
success_flag = false;

%%  PRIORITY 1: LOOKUP IN CUSTOM COORDINATES FILE
if ~isempty(custom_file) && exist(custom_file, 'file')
    [~, rinex_filename, ~] = fileparts(input_files.rinex_obs);
    if numel(rinex_filename) < 4
        fprintf('Warning: RINEX filename "%s" is too short to extract a 4-char station ID.\n', rinex_filename);
    else
        station_name_to_find = upper(rinex_filename(1:4));
        [found_xyz] = parse_custom_coords(custom_file, station_name_to_find);

        if ~isempty(found_xyz)
            final_xyz = found_xyz;
            success_flag = true;
            return; 
        else
        end
    end
end


%  PRIORITY 2: gLAB PPP CALCULATION 
if ~exist(glab_paths.executable, 'file')
    message = sprintf('gLAB executable not found at: %s. Fallback to approx_xyz.', glab_paths.executable);
    return;
end
if ~exist(glab_paths.working_dir, 'dir')
    [status_mkdir, msg_mkdir] = mkdir(glab_paths.working_dir);
    if ~status_mkdir
        message = sprintf('Could not create working directory: %s. Reason: %s. Fallback to approx_xyz.', glab_paths.working_dir, msg_mkdir);
        return;
    end
end

% --- Build gLAB command ---
base_cmd = sprintf('cd "%s" && "%s"', glab_paths.working_dir, glab_paths.executable);
cmd_parts = {base_cmd};
cmd_parts{end+1} = sprintf('-input:obs "%s"', input_files.rinex_obs);
if isfield(input_files, 'sp3') && ~isempty(input_files.sp3), cmd_parts{end+1} = sprintf('-input:sp3 "%s"', input_files.sp3); end
if isfield(input_files, 'ant') && ~isempty(input_files.ant), cmd_parts{end+1} = sprintf('-input:ant "%s"', input_files.ant); end
common_params = {'-pre:dec 300 -print:output', '-pre:elevation 30', '-model:solidtides', '-filter:nav static', '-filter:trop', '--print:info --print:cycleslips --print:input --print:model --print:satellites --print:prefit --print:postfit --print:filter --print:satdiff --print:satstat --print:satstattot --print:satpvt --print:sbasout --print:usererror --print:dgnss --print:summary --print:progress'};
cmd_parts = [cmd_parts, common_params];
command_str_attempt1 = strjoin(cmd_parts, ' ');

[~, cmdout1] = system(command_str_attempt1);

last_valid_xyz = parse_glab_output(cmdout1);

if ~isempty(last_valid_xyz)
    final_xyz = last_valid_xyz;
    success_flag = true;
    return;

elseif contains(cmdout1, "ERROR Reference station antenna")
    cmd_parts_retry = [cmd_parts, {'--model:recphasecenter'}];
    command_str_attempt2 = strjoin(cmd_parts_retry, ' ');
    [~, cmdout2] = system(command_str_attempt2);

    last_valid_xyz_retry = parse_glab_output(cmdout2);

    if ~isempty(last_valid_xyz_retry)
        final_xyz = last_valid_xyz_retry;
        success_flag = true;
        return; 
    else
        return;
    end
else

    return;
end
end

%  HELPER FUNCTION TO PARSE CUSTOM COORDINATE FILE
function [xyz] = parse_custom_coords(file_path, station_to_find)
xyz = [];
    fid = fopen(file_path, 'r');
    if fid == -1; return; end
    raw_text = fscanf(fid, '%c');
    fclose(fid);

    clean_text = regexprep(strrep(raw_text, newline, ' '), '\s+', ' ');
    parts = strsplit(strtrim(clean_text), ' ');

    for i = 1:numel(parts)
        if strcmpi(parts{i}, station_to_find)
            if i + 3 <= numel(parts)
                x = str2double(parts{i+1});
                y = str2double(parts{i+2});
                z = str2double(parts{i+3});
                if ~isnan(x) && ~isnan(y) && ~isnan(z)
                    xyz = [x, y, z];
                    return;
                end
            end
        end
    end

end


% --- Helper function to parse coordinates from gLAB's console output ---
function xyz = parse_glab_output(cmdout)
xyz = []; 
lines = splitlines(cmdout);
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if startsWith(line, 'OUTPUT')
        parts = strsplit(line);
        if numel(parts) >= 13
            x = str2double(parts{12});
            y = str2double(parts{13});
            z = str2double(parts{14});
            if ~isnan(x) && ~isnan(y) && ~isnan(z)
                xyz = [x, y, z]; 
            end
        end
    end
end
end

% Check and remove splice comments and last blank lines from RINEX files
function removeRinexSplice(filePath)

if ~isstring(filePath) && ~ischar(filePath)
    error('File path must be a string or character vector.');
end

if ~isfile(filePath)
    error('File not found: %s', filePath);
end


%% Define target string and read the file
targetString = 'RINEX FILE SPLICE; other post-header comments skipped       COMMENT';
lines = readlines(filePath);
if isempty(lines)
    return;
end

%% Find and mark lines for deletion
linesToDelete = false(size(lines));
spliceFound = false;
trailingBlanksFound = false;

for i = 2:numel(lines)
    if contains(lines(i), targetString)
        linesToDelete(i) = true;
        linesToDelete(i-1) = true;
        spliceFound = true;
    end
end

lastIndex = numel(lines);
while lastIndex > 0 && strtrim(lines(lastIndex)) == ""
    if ~linesToDelete(lastIndex)
        linesToDelete(lastIndex) = true;
        trailingBlanksFound = true;
    end
    lastIndex = lastIndex - 1;
end

if ~spliceFound && ~trailingBlanksFound
    return;
end

keptLines = lines(~linesToDelete);


try
    fileID = fopen(filePath, 'w');
    cleanupObj = onCleanup(@() fclose(fileID));

    if ~isempty(keptLines)
        fprintf(fileID, '%s\n', keptLines);
    end
      
catch ME
    error('Error writing to file: %s\nError message: %s', filePath, ME.message);
end

end
