function [] = simulation_scan_4wave(min,max)
    for l = min:max
        PIA_Gain(l,100,1310)
    end
end
