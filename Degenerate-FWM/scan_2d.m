function [] = scan_2d(min,max,min_pump,max_pump)
    for p = min_pump:5:max_pump
        simulation_scan_4wave(min,max,p)
    end
end