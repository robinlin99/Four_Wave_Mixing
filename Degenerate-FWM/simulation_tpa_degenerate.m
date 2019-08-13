function [] = simulation_tpa_degenerate(l,pl,power)
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
    [x_vec1,s_gain1,sig_gain1] = PIA_Gain(l,power,pl,0);
    [x_vec2,s_gain2,sig_gain2] = PIA_Gain(l,power,pl,1);
    [x_vec3,s_gain3,sig_gain3] = PIA_Gain(l,power,pl,2);
    [x_vec4,s_gain4,sig_gain4] = PIA_Gain(l,power,pl,3);

    % plot for different values of tpa 
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
    fname = '/Users/robinlin/Desktop/Research/2019 Summer/Four_Wave_Mixing/Experiment_scanner_tpa_degenerate';
    fng = sprintf('%0.2f nm Signal Gain degenerate tpa.png',l);
    saveas(gcf, fullfile(fname,fng));
end
