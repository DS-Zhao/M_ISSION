function [freq,wavl] = frequencies
%% Calculate Frequencies and Wavelengths for GNSS Satellites
% Function:
%     Generates carrier frequencies and corresponding wavelengths for various 
%     Global Navigation Satellite Systems (GNSS). This function contains 
%     hard-coded nominal values for the primary frequency bands 
%     (e.g., L1/E1/B1 and L2/E5a/B2) of GPS, GLONASS, Galileo, and BeiDou satellites.
%     It takes no input arguments and returns two matrices containing the 
%     frequencies (in Hz) and wavelengths (in meters).
%
% INPUT:
%     None
%
% OUTPUT:
%     freq:       An N x 2 matrix containing the satellite frequencies in Hertz (Hz). 
%                 Column 1 typically corresponds to the L1/E1/B1 band, and 
%                 Column 2 to the L2/E5a/B2 band.
%     wavl:       An N x 2 matrix containing the corresponding wavelengths in meters (m).
%% ---------------------------------------------------------------------

c = 299792458;
freq = zeros(92,2);
wavl = zeros(92,2);
glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];
for i=1:105
    if i<33     % GPS 32
        freq(i,1) = 10.23*10^6*154;    %Hz
        wavl(i,1) = c/(10.23*10^6*154);%m
        freq(i,2) = 10.23*10^6*120;    %Hz
        wavl(i,2) = c/(10.23*10^6*120);%m
    elseif i<59 % GLONASS 26
        freq(i,1) = (1602 + 0.5625*glok(i-32))*10^6;    %Hz
        wavl(i,1) = c/((1602 + 0.5625*glok(i-32))*10^6);%m
        freq(i,2) = (1246 + 0.4375*glok(i-32))*10^6;    %Hz
        wavl(i,2) = c/((1246 + 0.4375*glok(i-32))*10^6);%m
    elseif i<89 % GALILEO 30
        freq(i,1) = 10.23*10^6*154;    %Hz
        wavl(i,1) = c/(10.23*10^6*154);%m
        freq(i,2) = 10.23*10^6*115;    %Hz
        wavl(i,2) = c/(10.23*10^6*115);%m
    else        % BEIDOU 16
        freq(i,1) = 10.23*10^6*152.6;    %Hz
        wavl(i,1) = c/(10.23*10^6*152.6);%m
        freq(i,2) = 10.23*10^6*118;
        wavl(i,2) = c/(10.23*10^6*118);%m
    end
end
end

