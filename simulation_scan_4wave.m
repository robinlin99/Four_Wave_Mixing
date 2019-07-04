function [] = simulation_scan_4wave(min,max,pump)
    % Initialize
    sig_gain_arr = zeros(1,max-min+1);
    index = 1;
    sig_gain = 0;
    wavelength = linspace(min,max,max-min+1);
    for l = min:max
        [sig_gain] = PIA_Gain(l,pump,1310);
        sig_gain_arr(index) = sig_gain;
        index = index + 1;
    end
    figure(1)
    plot(wavelength,sig_gain_arr);
    xlabel 'Wavelength (nm)'
    ylabel 'Gain (dB)'
    grid on
    grid minor
    fname = '/Users/robinlin/Desktop/Research/2019 Summer/Four_Wave_Mixing/Experiment_Scanner'
    fng = sprintf('Pump = %0.2f mW Signal Gain.eps', pump);
    saveas(gcf, fullfile(fname,fng));
end
