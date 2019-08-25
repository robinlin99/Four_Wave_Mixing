function [] = simulation_scan_4wave_4wave_spm(min,max)
    %% Initialize
    % We will do the scan at four different pump power in order to
    % visualize the effects of self-phase modulation 
    % pump power at 30,40,50,60
    sig_gain_arr30 = zeros(1,max-min+1);
    sig_gain_arr40 = zeros(1,max-min+1);
    sig_gain_arr50 = zeros(1,max-min+1);
    sig_gain_arr60 = zeros(1,max-min+1);
    index = 1;
    sig_gain30 = 0;
    sig_gain40 = 0;
    sig_gain50 = 0;
    sig_gain60 = 0;
    wavelength = linspace(min,max,max-min+1);
    x_vec30 = zeros(1,101);
    x_vec40 = zeros(1,101);
    x_vec50 = zeros(1,101);
    x_vec60 = zeros(1,101);
    s_gain30 = zeros(1,101);
    s_gain40 = zeros(1,101);
    s_gain50 = zeros(1,101);
    s_gain60 = zeros(1,101);
    for l = min:max
        [x_vec30,s_gain30,sig_gain30] = PIA_Gain_4wave(l,30,30,1310,1310,0);
        [x_vec40,s_gain40,sig_gain40] = PIA_Gain_4wave(l,40,40,1310,1310,0);
        [x_vec50,s_gain50,sig_gain50] = PIA_Gain_4wave(l,50,50,1310,1310,0);
        [x_vec60,s_gain60,sig_gain60] = PIA_Gain_4wave(l,60,60,1310,1310,0);
        sig_gain_arr30(index) = sig_gain30;
        sig_gain_arr40(index) = sig_gain40;
        sig_gain_arr50(index) = sig_gain50;
        sig_gain_arr60(index) = sig_gain60;
        index = index + 1;
    end
    figure(1)
    plot(wavelength,sig_gain_arr30);
    hold on
    plot(wavelength,sig_gain_arr40);
    hold on
    plot(wavelength,sig_gain_arr50);
    hold on
    plot(wavelength,sig_gain_arr60);
    xlabel 'Wavelength (nm)'
    ylabel 'Gain (dB)'
    legend('30 mW','40 mW', '50 mW','60 mW');
    grid on
    grid minor
    fname = '/Users/robinlin/Desktop/Research/2019 Summer/Four_Wave_Mixing/Experiment_Scanner_Nondegenerate_spm';
    fng = sprintf('SPM Nondegenerate Signal Gain.png');
    saveas(gcf, fullfile(fname,fng));
end
