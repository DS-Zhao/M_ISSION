function [data] = preprocess(data,options)

[data]  = elm_badclk(data);%Clear observation data with satellite coordinate errors

[data] = decimation(data,options);%Clear observation data that is not within the calculation range

[data] = cal_sat(data);%Calculate the actual satellite coordinates

[data] = elv_mask(data,options);%Cut off angle

[data] = cs_detect(data,options);%Cycle jump detection and repair, optional GF method or MW method

if options.clkjump == 1
    [data] = clk_jmp2(data);%Receiver clock jump
end

[data] = outlier(data);%Exclude observation data with gross errors

if options.codsmth == 1
    [data] = smoothing(data);%Carrier pseudorange smoothing
end
end

