function [] = simulation_tpa(l,pump1_l,pump2_l,pump1,pump2)
    sig_gain1 = 0;
    sig_gain2 = 0;
    sig_gain3 = 0;
    sig_gain4 = 0;
    x_vec1 = zeros(1,101);
    x_vec2 = zeros(1,101);
    x_vec3 = zeros(1,101);
    x_vec4 = zeros(1,101);
    s_gain1 = zeros(1,101);
    s_gain2 = zeros(1,101);
    s_gain3 = zeros(1,101);
    s_gain4 = zeros(1,101);
    [x_vec1,s_gain1,sig_gain1] = PIA_Gain_4wave(l,pump1, pump2,pump1_l,pump2_l,0);
    [x_vec2,s_gain2,sig_gain2] = PIA_Gain_4wave(l,pump1, pump2,pump1_l,pump2_l,1);
    [x_vec3,s_gain3,sig_gain3] = PIA_Gain_4wave(l,pump1, pump2,pump1_l,pump2_l,2);
    [x_vec4,s_gain4,sig_gain4] = PIA_Gain_4wave(l,pump1, pump2,pump1_l,pump2_l,3);

    % plot no two photon absorption with two photon absorption on
    figure;
    plot(x_vec1,s_gain1);
    hold on
    plot(x_vec2,s_gain2);
    hold on
    plot(x_vec3,s_gain3);
    hold on
    plot(x_vec4,s_gain4);
    grid on
    grid minor
    xlabel 'Length (mm)';
    ylabel 'Gain (dB)';
    legend('TPA - 1.5e-12','TPA - 2e-12', 'TPA - 2.5e-12','TPA - 3e-12');
    fname = '/Users/robinlin/Desktop/Research/2019 Summer/Four_Wave_Mixing/Experiment_scanner_tpa';
    fng = sprintf('%0.2f nm, Pump1 = %0.2f mW, Pump2 = %0.2f mW Signal Gain.png',l, pump1, pump2);
    saveas(gcf, fullfile(fname,fng));
end
