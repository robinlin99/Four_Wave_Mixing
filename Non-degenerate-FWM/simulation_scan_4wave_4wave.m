function [] = simulation_scan_4wave_4wave(min,max,pump1_l,pump2_l,pump1,pump2)
    for l = min:max
        PIA_Gain_4wave(l, pump1, pump2, pump1_l, pump2_l);
    end
end
