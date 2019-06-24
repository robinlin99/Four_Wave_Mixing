function [] = simulation_scan_4wave_nd(min,max)
% Scans from minimum to maximum wavelength at a step length of 1nm 
% Keeping pump wavelength at 1310 nm
% Need to find the phase-matching condition
% No output
    for l = min:max
        PIA_Gain_nd(l,100,1310)
    end
    
    
end

