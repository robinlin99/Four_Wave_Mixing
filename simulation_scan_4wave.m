function [] = simulation_scan_4wave(min,max)
% Scans from minimum to maximum wavelength at a step length of 1nm 
% No output
    for l = min:max
        PIA_Gain(l,20,1310)
    end
end

