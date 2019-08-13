function [] = simulation_spm_degenerate(l,pump)
    sig_gain1 = 0;
    sig_gain2 = 0;
    sig_gain3 = 0;
    sig_gain4 = 0;
    sig_gain5 = 0;
    sig_gain6 = 0;
    x_vec1 = zeros(1,101);
    x_vec2 = zeros(1,101);
    x_vec3 = zeros(1,101);
    x_vec4 = zeros(1,101);
    x_vec5 = zeros(1,101);
    x_vec6 = zeros(1,101);
    s_gain1 = zeros(1,101);
    s_gain2 = zeros(1,101);
    s_gain3 = zeros(1,101);
    s_gain4 = zeros(1,101);
    s_gain5 = zeros(1,101);
    s_gain6 = zeros(1,101);
    [x_vec1,s_gain1,sig_gain1] = PIA_Gain(l,30,pump,0);
    [x_vec2,s_gain2,sig_gain2] = PIA_Gain(l,40,pump,0);
    [x_vec3,s_gain3,sig_gain3] = PIA_Gain(l,50,pump,0);
    [x_vec4,s_gain4,sig_gain4] = PIA_Gain(l,60,pump,0);
    [x_vec5,s_gain5,sig_gain5] = PIA_Gain(l,70,pump,0);
    [x_vec6,s_gain6,sig_gain6] = PIA_Gain(l,80,pump,0);

    % plot for different values of tpa 
    figure;
    plot(x_vec1,s_gain1);
    hold on
    plot(x_vec2,s_gain2);
    hold on
    plot(x_vec3,s_gain3);
    hold on
    plot(x_vec4,s_gain4);
    hold on
    plot(x_vec5,s_gain5);
    hold on
    plot(x_vec6,s_gain6);
    grid on
    grid minor
    xlabel 'Length (mm)';
    ylabel 'Gain (dB)';
    legend('Power - 30 mW','Power - 40 mW', 'Power - 50 mW','Power - 60 mW',"Power - 70 mW","Power - 80 mW");
    fname = '/Users/robinlin/Desktop/Research/2019 Summer/Four_Wave_Mixing/Experiment_scanner_spm_degenerate';
    fng = sprintf('%0.2f nm Signal Gain degenerate spm.png',l);
    saveas(gcf, fullfile(fname,fng));
end
