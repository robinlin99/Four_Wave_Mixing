function rst = PIA_Gain_nd(wl, pump, pl)

%  2-D BPM code
%  chi3 (quadratic) nonlinear material
%  PIA_Gain(wl, pump, pl, isTE)

format long
simm=1;
t0=clock;
j=sqrt(-1);   % imaginary unit
PSA = 1;
%wlstep=6;      % signal and idler wavelength difference in the unit of nm

%% Assignment for Idler, Signal, and Pump

% 1 - Signal
% 2 - Idler
% 3 - Pump
%% input data
c=299792458   ;                 % velocity of light
epsilon = 8.85e-12;
% Symbolic library
syms x w
% BRW waveguide TE and TM mode refractive index [x is wavelength in meter]
% Relevant parameters derived from the refractive index profile
N_te(x) = 0.1604e1 + 0.1107e7*x - 0.4101e12*x^2 + 0.5822e17*x^3 + 0.9198e-6/x;
N_tm(x) = 1.28 + 0.1416e7*x - 0.5435e12*x^2 + 0.7951e17 * x^3 + 0.1039e-5/x;
N_te_w(w) = subs(N_te,x,2*pi*c/w);
N_tm_w(w) = subs(N_tm,x,2*pi*c/w);
d_te(w) = diff(N_te_w,w);
d_tm(w) = diff(N_tm_w,w);
b_1_te(w) = (1/c)*(N_te_w+w*d_te);
b_1_te_l(x) = subs(b_1_te,w,2*pi*c/x);
b_2_te(w) = (1/c)*(2*d_te+w*diff(d_te));
b_2_te_l(x) = subs(b_2_te,w,2*pi*c/x);
b_1_tm(w) = (1/c)*(N_tm_w+w*d_tm);
b_1_tm_l(x) = subs(b_1_tm,w,2*pi*c/x);
b_2_tm(w) = (1/c)*(2*d_tm+w*diff(d_tm));
b_2_tm_l(x) = subs(b_2_tm,w,2*pi*c/x);
v_g_te(w) = c/(N_te_w+w*d_te);
v_g_te_l(x) = subs(v_g_te,w,2*pi*c/x);
v_g_tm(w) = c/(N_tm_w+w*d_tm);
v_g_tm_l(x) = subs(v_g_tm,w,2*pi*c/x);
gvd_te(x) = subs(diff(1/v_g_te),w,2*pi*c/x);
gvd_tm(x) = subs(diff(1/v_g_tm),w,2*pi*c/x);
% Substitute wavelength back in to the equation
b_1_te_p(x) = subs(b_1_te,w,2*pi*c/x);
b_2_te_p(x) = subs(b_2_te,w,2*pi*c/x);
b_1_tm_p(x) = subs(b_1_tm,w,2*pi*c/x);
b_2_tm_p(x) = subs(b_2_tm,w,2*pi*c/x);

%-------------------------------------------------------
%TE AND TM
%pmp = load('SH3_TE_pump.mat');
%sig = load('SH3_TE_Signal_40nm.mat');
%idl = load('SH3_TM_idler1540-80nm.mat');

%----------------------------------------------------------------------
% Process | A_1 | A_2 | A_3 | A_4 | Frequency shift | Condition
%    I    |  s  |  f  |  s  |  f  |                 |   beta_2 > 0
%    II   |  s  |  f  |  f  |  s  |                 |   beta_2 < 0
%    III  |  s  |  s  |  f  |  f  |                 |   beta_2 > 0
%    IV   |  f  |  f  |  s  |  s  |                 |   beta_2 < 0
%----------------------------------------------------------------------


% Pump is in TE mode (slow axis)
pmp_n = N_te;
pmp_gvd = gvd_te;
pmp_b1 = b_1_te_l;
pmp_b2 = b_2_te_l;
pmp_vg = v_g_te_l;

% Signal and idler are in TM mode (fast axis)
sig_n = N_tm;      % make all in TM mode
idl_n = N_tm;
sig_gvd = gvd_tm;
idl_gvd = gvd_tm;
sig_b1 = b_1_tm_l;
idl_b1 = b_1_tm_l;
sig_b2 = b_2_tm_l;
idl_b2 = b_2_tm_l;
sig_vg = v_g_tm_l;
idl_vg = v_g_tm_l;


%-------OLD CODE -------------------
%sig.GVD =-(sig.D).* c./((sig.f).^2)/2/pi;
%sig.l = c./(sig.f);
%sig.neff = real(sig.neff);

%idl.GVD =-(idl.D).* c./((idl.f).^2)/2/pi;
%idl.l = c./(idl.f);
%idl.neff = real(idl.neff);

%pmp.GVD = -(pmp.D).* c./((pmp.f).^2)/2/pi;
%pmp.l = c./(pmp.f);
%pmp.neff = real(pmp.neff);
%------OLD CODE -------------------

L=1*1.0e-3;                         % length of the waveguide
lambda_wave1= wl*1e-9;        % 1st wave wavelength (signal)
lambda_wave3= pl*1e-9; % 3rd wave wavelength (pump)
lambda_wave2 = 1/(2/lambda_wave3-1/lambda_wave1); % 2nd wavelength, from energy conservation, (idler)

k0_1=2*pi/lambda_wave1;        % free-space propagation constant 1st wave
k0_2=2*pi/lambda_wave2;        % free-space propagation constant 2nd wave
k0_3=2*pi/lambda_wave3;        % free-space propagation constant 3rd wave

%% refractive index data

%n_wave1_d=interp1(sig.l, sig.neff, lambda_wave1);        % SH3
%n_wave2_d=interp1(idl.l, idl.neff, lambda_wave2);
%n_wave3_d=interp1(pmp.l, pmp.neff, lambda_wave3);
n_wave1_d = double(sig_n(lambda_wave1));
n_wave2_d = double(idl_n(lambda_wave2));
n_wave3_d = double(pmp_n(lambda_wave3));
%% refractive indices in material as grown

n_wave1_o=n_wave1_d;
n_wave2_o=n_wave2_d;
n_wave3_o=n_wave3_d;

%% Propagation constants in the different media

k_wave1_d=2*pi*n_wave1_d/lambda_wave1; % 1st wave propagation constant in disordered material
k_wave2_d=2*pi*n_wave2_d/lambda_wave2; % 2nd wave propagation constant in disordered material
k_wave3_d=2*pi*n_wave3_d/lambda_wave3; % 3rd wave propagation constant in disordered material


%% Group Velocities in the different media

%vg_wave1_d = interp1(sig.l, sig.vg, lambda_wave1);  % for SH3
%vg_wave2_d = interp1(idl.l, idl.vg, lambda_wave2);
%vg_wave3_d = interp1(pmp.l, pmp.vg, lambda_wave3);

vg_wave1_d = double(sig_vg(lambda_wave1));
vg_wave2_d = double(idl_vg(lambda_wave2));
vg_wave3_d = double(pmp_vg(lambda_wave3));

%% Group Velocity Mismatch Parameters

tao = 1; % pulse duration

gamma_wave1_d = 1/vg_wave1_d;
gamma_wave2_d = 1/vg_wave2_d;
gamma_wave3_d = 1/vg_wave3_d;


N_wave1_d = 1*(-gamma_wave3_d+gamma_wave1_d)/tao;
N_wave1_o = N_wave1_d;

N_wave2_d = (-gamma_wave3_d+gamma_wave2_d)/tao;
N_wave2_o = N_wave2_d;

N_wave3_d = 1*(-gamma_wave3_d+gamma_wave3_d)/tao;
N_wave3_o = N_wave3_d;

%% Group Velocity Dispersion Parameters
%beta_wave1_o = interp1(sig.l, sig.GVD, lambda_wave1);                      %SH3
%beta_wave2_o = interp1(idl.l, idl.GVD, lambda_wave2);
%beta_wave3_o = interp1(pmp.l, pmp.GVD, lambda_wave3);

beta_wave1_o = double(sig_gvd(lambda_wave1));
beta_wave2_o = double(idl_gvd(lambda_wave2));
beta_wave3_o = double(pmp_gvd(lambda_wave3));

%% Linear Losses

alpha_wave3_o = 200;
alpha_wave2_o = 200;
alpha_wave1_o = 200;

%% n2

n2_wave1_o = 1*100*1e-20;            %SH3
n2_wave2_o = 1*100*1e-20;
n2_wave3_o = 600*1e-20;
%% alpha2 two photon absorption
a2_wave1_o = 1.5*1e-12;                 % SH3
a2_wave2_o = 1.5*1e-12;
a2_wave3_o = 1.5*1e-12;

%% The pump and signal power amplitude
gradual = 1;
T0=0.064*1e-12 /1.665;    % T0=FWHM/1.665
T0_wave1=100*T0;
T0_wave2=100*T0;
T0_wave3=T0;
f0 = 80e6; % repetition rate of pump laser


amp_wave3_sqar=(1.09375*0.2*pump*0.75*10e2)*2^(gradual-1);  %peak power for pulsed pump electric field is square root
amp_wave1_sqar=1e3*1e-12;
amp_wave2_sqar=1e-20;
x = 0:0.01:1;

%% Phase-Mismatch and Coherent Lengths

% K1 and K2 -> Idler
% K3 -> Pump
deltak_d = -2*k_wave3_d + k_wave1_d + k_wave2_d;
deltak_o = deltak_d;

lc_d = pi/deltak_d; % coherent length in disordered material
lc_o = lc_d; % coherent length in material as grown

%%
%% Number of Periods used

num_period =100;

%% Numerical Simulation Data

pointt= 1024;            % number of points along t
iterations=5;         % number of iterations per material
numplots=100;           % the number of plots along z is numplots+1
njump=(iterations*2*num_period)/numplots;  % indicator used to know when to save the data

tmax=1*30e-12;             % the spatial window is from -tmax to tmax

z_period=(lc_d+lc_o);   % size of a single period


deltat=2*tmax/pointt;   % step along t


%% temporary structure resolution per material

%%%%% this section may seem redundant, although it can be used to
%%%%% introduce arbitrary number of points per medium rather than 100 each
T =  [];  % structure pattern
Ssign = 1;
for zstep = 1:2*num_period*iterations
    T(zstep) = 1*Ssign;
    if rem(zstep,iterations) == 0
        Ssign = -Ssign;
    end
end

changer = []; % array used to track the material step

for st = 2:length(T)-1;
    if T(st) ~= T(st+1)
        changer = [changer st];
    end
end

changer = [0 changer length(T)]; %position of each boundary

%% GVD and GVM vectors
indfreq=-pointt/2:1:pointt/2-1; % number of point along t
omega=(pi./tmax).*indfreq;    % temporal angular frequency

%
gvd_wave1_o=1.0*(beta_wave1_o/2)*(omega.^2); % 1st wave group velocity dispersion
gvd_wave2_o=1.0*(beta_wave2_o/2)*(omega.^2); % 2nd wave group velocity dispersion
gvd_wave3_o=1.0*(beta_wave3_o/2)*(omega.^2); % 3rd wave group velocity dispersion

gvd_wave1_d=gvd_wave1_o;                 % 1st wave group velocity dispersion
gvd_wave2_d=gvd_wave2_o;                 % 2nd wave group velocity dispersion
gvd_wave3_d=gvd_wave3_o;                 % 3rd wave group velocity dispersion

gvm_wave1_d= 1.0*N_wave1_d*(omega);        % 1st wave group velocity mismatch
gvm_wave2_d= 1.0*N_wave2_d*(omega);        % 2nd wave group velocity mismatch
gvm_wave3_d= 1.0*N_wave3_d*(omega);        % 3rd wave group velocity mismatch

gvm_wave1_o= gvm_wave1_d;                % 1st wave group velocity mismatch
gvm_wave2_o= gvm_wave2_d;                % 2nd wave group velocity mismatch
gvm_wave3_o= gvm_wave3_d;                % 3rd wave group velocity  mismatch

gvm_wave1_m = [];
gvm_wave2_m = [];
gvm_wave3_m = [];

gvd_wave1_m = [];
gvd_wave2_m = [];
gvd_wave3_m = [];

%% Building the structure
chi2_1 = 70*1e-12 ;                      % (effective) quadratic nonlinear coefficient
chi2_2 = chi2_1;
A_eff_o = 1e-12;                   % SH3
A_eff_d =  A_eff_o;

A_eff_3_TE = 6.553*1e-12;               % SH3
A_eff_3_TM = 6.665*1e-12;
A_eff_3_o = A_eff_o;

%A_eff_3_d = A_eff_3_o;

Tr_wave3(1:num_period*2*iterations+2) = 1;
Tr_wave2(1:num_period*2*iterations+2) = 1;
Tr_wave1(1:num_period*2*iterations+2) = 1;

if simm == 1
%     lc_1 = lc_o;
      lc_1 = L/(2*num_period);
      lc_2 = lc_1;
elseif simm == 2
    lc_1 = 0.5*z_period;
    lc_2 = 0.5*z_period;
elseif simm == 3
    lc_1 = 0.35*z_period;
    lc_2 = 0.65*z_period;
elseif simm ==4
    lc_1 = K_lc1*z_period;
    lc_2 = (1-K_lc1)*z_period;
end

zStepL = length(changer); % this is how many steps in the z-direction for differential equation solver
z_m(1:zStepL-1) = [0];
n_wave1_m(1:zStepL-1) = [0];
n_wave2_m(1:zStepL-1) = [0];
n_wave3_m(1:zStepL-1) = [0];
deltak_m(1:zStepL-1) = [0];
chi2_m(1:zStepL-1) = [0];
gvm_wave1_m(1:zStepL-1,length(gvm_wave1_d)) = [0];
gvm_wave2_m(1:zStepL-1,length(gvm_wave2_d)) = [0];
gvm_wave3_m(1:zStepL-1,length(gvm_wave3_d)) = [0];
alpha_wave1_m(1:zStepL-1) = [0];
alpha_wave2_m(1:zStepL-1) = [0];
alpha_wave3_m(1:zStepL-1) = [0];
A_eff_m(1:zStepL-1) = [0];
n2_wave1_m(1:zStepL-1) = [0];
n2_wave2_m(1:zStepL-1) = [0];
n2_wave3_m(1:zStepL-1) = [0];
a2_wave1_m(1:zStepL-1) = [0];
a2_wave2_m(1:zStepL-1) = [0];
a2_wave3_m(1:zStepL-1) = [0];
gvd_wave1_m(1:zStepL-1,length(gvd_wave1_d)) = [0];
gvd_wave2_m(1:zStepL-1,length(gvd_wave2_d)) = [0];
gvd_wave3_m(1:zStepL-1,length(gvd_wave3_d)) = [0];
A_eff_3_m(1:zStepL-1) = [0];

for zstep = 1:zStepL-1
    z_m(zstep) = lc_1;
    n_wave1_m(zstep) = n_wave1_o;
    n_wave2_m(zstep) = n_wave2_o;
    n_wave3_m(zstep) = n_wave3_o;

    deltak_m(zstep) = deltak_o;
    chi2_m(zstep) = chi2_1;

    gvm_wave1_m(zstep,:) = gvm_wave1_o;
    gvm_wave2_m(zstep,:) = gvm_wave2_o;
    gvm_wave3_m(zstep,:) = gvm_wave3_o;

    gvd_wave1_m(zstep,:) = gvd_wave1_o;
    gvd_wave2_m(zstep,:) = gvd_wave2_o;
    gvd_wave3_m(zstep,:) = gvd_wave3_o;

    alpha_wave3_m(zstep) = alpha_wave3_o;
    alpha_wave2_m(zstep) = alpha_wave2_o;
    alpha_wave1_m(zstep) = alpha_wave1_o;

    A_eff_m(zstep) = A_eff_o;

    n2_wave1_m(zstep) = n2_wave1_o;
    n2_wave2_m(zstep) = n2_wave2_o;
    n2_wave3_m(zstep) = n2_wave3_o;

    a2_wave1_m(zstep) = a2_wave1_o;
    a2_wave2_m(zstep) = a2_wave2_o;
    a2_wave3_m(zstep) = a2_wave3_o;

    A_eff_3_m(zstep) = A_eff_3_o;
end

length_m = changer(2:end) - changer(1:end-1);

deltaz_m = z_m./5;
% The complete Structure with the different properties for each material

Structure = {length_m deltaz_m n_wave1_m n_wave2_m n_wave3_m deltak_m chi2_m gvm_wave1_m gvm_wave2_m gvm_wave3_m gvd_wave1_m gvd_wave2_m gvd_wave3_m alpha_wave3_m alpha_wave2_m alpha_wave1_m A_eff_m n2_wave3_m n2_wave2_m n2_wave1_m a2_wave1_m a2_wave2_m a2_wave3_m A_eff_3_m};

%%  Input Field
amp_wave1=0;
amp_wave2=0;
amp_wave3=0;
% amp_wave1_sqar=0;
% amp_wave2_sqar=0;
% amp_wave3_sqar=0;
%% Loop for input power
gradual_max=1;

power_wave1_L=[];
power_wave2_L=[];
power_wave3_L=[];
eta_L=[];

for gradual=1:1:gradual_max

amp_wave1=sqrt(amp_wave1_sqar);   % input 1st wave amplitude
amp_wave2=sqrt(amp_wave2_sqar);   % input 2nd wave amplitude
amp_wave3=sqrt(amp_wave3_sqar);   % input 3rd wave amplitude

t = -tmax:deltat:tmax-deltat;     % time interval

%% Pulse shape and initilizations, Gaussian shape is generally used, but sech() wave is better for femtosecond pulses

e_wave1=amp_wave1*exp(-(t/(sqrt(2)*T0_wave1)).^2);   % input 1st wave Gaussian beam
e_wave2=amp_wave2*exp(-(t/(sqrt(2)*T0_wave2)).^2);   % input 2nd wave Gaussian beam
e_wave3=amp_wave3*sech(-(t/(sqrt(1)*T0_wave3)).^1);   % ultrafast pulse in sech form
%e_wave3=amp_wave3*exp(-(t/(sqrt(2)*T0_wave3)).^2);   % input 3rd wave Gaussian beam

%% phase and z direction starting point and initial values
zeta = 0; %z starting point
pphase = 0; %phase difference starting point

%% vector initialization

distance=[];
i_wave1=[];           % 1st wave intensity
i_wave2=[];           % 2nd wave intensity
i_wave3=[];           % 3rd wave intensity
inorm_wave1=[];       % 1st wave normalized intensity
inorm_wave2=[];       % 2nd wave normalized intensity
inorm_wave3=[];       % 3rd wave normalized intensity
e_wave1_field=[];     % 1st wave field
e_wave2_field=[];     % 2nd wave field
e_wave3_field=[];     % 3rd wave field

i_fft_wave1=[];       % 1st wave power spectrum
i_fft_wave2=[];       % 2nd wave power spectrum
i_fft_wave3=[];       % 3rd wave power spectrum

i_fftnorm_wave1=[];   % 1st wave normalized power spectrum
i_fftnorm_wave2=[];   % 2nd wave normalized power spectrum
i_fftnorm_wave3=[];   % 3rd wave normalized power spectrum

distance(1:numplots+1)=[0];
phaseplot= [];

%% Memory Pre-Allocation
      power_wave1(1:numplots+1)=[deltat*sum((abs(e_wave1).^2))] ;
      i_wave1(1:length(t),1:numplots+1)=[0];
      i_wave1(1:length(t),1)=[abs(e_wave1').^2];
      inorm_wave1(1:length(t),1:numplots+1)=[0];
      inorm_wave1(1:length(t),1)=[(abs(e_wave1').^2)/max(abs(e_wave1').^2)];
      e_wave1_field(1:length(t),1:numplots+1)=[0];
      e_wave1_field(1:length(t),1)=[e_wave1'];
      i_fft_wave1(1:length(t),1:numplots+1)=[0];
      i_fft_wave1(1:length(t),1)=[abs(deltat*fftshift(fft(e_wave1'))).^2 ];
      i_fftnorm_wave1(1:length(t),1:numplots+1)=[0];
      i_fftnorm_wave1(1:length(t),1)=[abs(fftshift(fft(e_wave1'))).^2/max(abs(fftshift(fft(e_wave1'))).^2) ];

      power_wave2(1:numplots+1)=[deltat*sum((abs(e_wave2).^2))] ;
      i_wave2(1:length(t),1:numplots+1)=[0];
      i_wave2(1:length(t),1)=[abs(e_wave2').^2];
      inorm_wave2(1:length(t),1:numplots+1)=[0];
      inorm_wave2(1:length(t),1)=[(abs(e_wave2').^2)/max(abs(e_wave2').^2)];
      e_wave2_field(1:length(t),1:numplots+1)=[0];
      e_wave2_field(1:length(t),1)=[e_wave2'];
      i_fft_wave2(1:length(t),1:numplots+1)=[0];
      i_fft_wave2(1:length(t),1)=[abs(deltat*fftshift(fft(e_wave2'))).^2 ];
      i_fftnorm_wave2(1:length(t),1:numplots+1)=[0];
      i_fftnorm_wave2(1:length(t),1)=[abs(fftshift(fft(e_wave2'))).^2/max(abs(fftshift(fft(e_wave2'))).^2) ];

      power_wave3(1:numplots+1)=[deltat*sum((abs(e_wave3).^2))] ;
      i_wave3(1:length(t),1:numplots+1)=[0];
      i_wave3(1:length(t),1)=[abs(e_wave3').^2];
      inorm_wave3(1:length(t),1:numplots+1)=[0];
      inorm_wave3(1:length(t),1)=[(abs(e_wave3').^2)/max(abs(e_wave3').^2)];
      e_wave3_field(1:length(t),1:numplots+1)=[0];
      e_wave3_field(1:length(t),1)=[e_wave3'];
      i_fft_wave3(1:length(t),1:numplots+1)=[0];
      i_fft_wave3(1:length(t),1)=[abs(deltat*fftshift(fft(e_wave3'))).^2 ];
      i_fftnorm_wave3(1:length(t),1:numplots+1)=[0];
      i_fftnorm_wave3(1:length(t),1)=[abs(fftshift(fft(e_wave3'))).^2/max(abs(fftshift(fft(e_wave3'))).^2) ];
%% MAIN PROGRAM
z_count = 0;
Material_step = 2; %material number indicator

for loop_step=1:1:num_period*2*iterations

deltaz = Structure{2}(Material_step - 1);
n_wave1 = Structure{3}(Material_step - 1);
n_wave2 = Structure{4}(Material_step - 1);
n_wave3 = Structure{5}(Material_step - 1);
deltak = Structure{6}(Material_step - 1);
chi2 = Structure{7}(Material_step - 1);
gvm_wave1 = Structure{8}((Material_step - 1),:);
gvm_wave2 = Structure{9}((Material_step - 1),:);
gvm_wave3 = Structure{10}((Material_step - 1),:);
gvd_wave1 = Structure{11}((Material_step - 1),:);
gvd_wave2 = Structure{12}((Material_step - 1),:);
gvd_wave3 = Structure{13}((Material_step - 1),:);
alpha_wave3 = Structure{14}(Material_step - 1);
alpha_wave2 = Structure{15}(Material_step - 1);
alpha_wave1 = Structure{16}(Material_step - 1);
A_eff_2nd = Structure{17}(Material_step - 1);
n2_wave3 = Structure{18}(Material_step - 1);
n2_wave2 = Structure{19}(Material_step - 1);
n2_wave1 = Structure{20}(Material_step - 1);
a2_wave1 = Structure{21}(Material_step - 1);
a2_wave2 = Structure{22}(Material_step - 1);
a2_wave3 = Structure{23}(Material_step - 1);
A_eff_3rd = Structure{24}(Material_step - 1);

% coupling coefficients
%---------------------------------------------------------------------------------------------------
% CHI2
%c_chi2_wave1=-j*chi2*sqrt(8*pi^2/(n_wave1*n_wave2*n_wave3*c*epsilon*lambda_wave1^2*A_eff_2nd));
%c_chi2_wave2=-j*chi2*sqrt(8*pi^2/(n_wave1*n_wave2*n_wave3*c*epsilon*lambda_wave2^2*A_eff_2nd));
%c_chi2_wave3=-j*chi2*sqrt(8*pi^2/(n_wave1*n_wave2*n_wave3*c*epsilon*lambda_wave3^2*A_eff_2nd));

%CHI3
c_chi3_wave1 = 2*j*(n2_wave1*2*pi)/(A_eff_3rd*lambda_wave1);
c_chi3_wave2 = 2*j*(n2_wave2*2*pi)/(A_eff_3rd*lambda_wave2);
c_chi3_wave3 = 2*j*(n2_wave3*2*pi)/(A_eff_3rd*lambda_wave3);

c_n2a2_wave1 = (j*2*pi*n2_wave1/lambda_wave1 - a2_wave1/2)/(A_eff_3rd);
c_n2a2_wave2 = (j*2*pi*n2_wave2/lambda_wave2 - a2_wave2/2)/(A_eff_3rd);
c_n2a2_wave3 = (j*2*pi*n2_wave3/lambda_wave3 - a2_wave3/2)/(A_eff_3rd);


% transmission coefficients
e_wave1 = e_wave1*Tr_wave1(loop_step);
e_wave2 = e_wave2*Tr_wave2(loop_step);
e_wave3 = e_wave3*Tr_wave3(loop_step);


% function that calculates the values of the pulses after each step
[e_wave1,e_wave2,e_wave3,es_wave1,es_wave2,es_wave3] = step_4wave(e_wave1,e_wave2,e_wave3,deltat,deltaz,deltak,gvd_wave1,gvd_wave2,gvd_wave3,c_chi3_wave1,c_chi3_wave2,c_chi3_wave3,gvm_wave1,gvm_wave2,gvm_wave3,pphase,alpha_wave1,alpha_wave2,alpha_wave3,c_n2a2_wave1,c_n2a2_wave2,c_n2a2_wave3);

pphase = deltak*deltaz + pphase; %% increment in the phase angle

% array of phase increments, this is also used to track the material steps
if loop_step == changer(Material_step)
    phaseplot = [phaseplot pphase];
    Material_step = 1 + Material_step;
end

zeta = deltaz + zeta; %increment in z

% Saving the results

   if rem(loop_step,njump)==0  % it saves the result (only if loop_step is multiple of njump)
      z_count = z_count+1;

      power_wave1(z_count+1)=[deltat*sum((abs(e_wave1).^2))] ;
      i_wave1(1:length(t),z_count+1)=[abs(e_wave1').^2];
      inorm_wave1(1:length(t),z_count+1)=[(abs(e_wave1').^2)/max(abs(e_wave1').^2)];
      e_wave1_field(1:length(t),z_count+1)=[e_wave1'];
      i_fft_wave1(1:length(t),z_count+1)=[abs(es_wave1').^2];
      i_fftnorm_wave1(1:length(t),z_count+1)=[abs(es_wave1').^2/max(abs(es_wave1').^2)];


      %
      power_wave2(z_count+1)=[deltat*sum((abs(e_wave2).^2))] ;
      i_wave2(1:length(t),z_count+1)=[abs(e_wave2').^2];
      inorm_wave2(1:length(t),z_count+1)=[(abs(e_wave2').^2)/max(abs(e_wave2').^2)];
      e_wave2_field(1:length(t),z_count+1)=[e_wave2'];
      i_fft_wave2(1:length(t),z_count+1)=[abs(es_wave2').^2];
      i_fftnorm_wave2(1:length(t),z_count+1)=[abs(es_wave2').^2/max(abs(es_wave2').^2)];


      %
      power_wave3(z_count+1)=[deltat*sum((abs(e_wave3).^2))] ;
      i_wave3(1:length(t),z_count+1)=[abs(e_wave3').^2];
      inorm_wave3(1:length(t),z_count+1)=[(abs(e_wave3').^2)/max(abs(e_wave3').^2)];
      e_wave3_field(1:length(t),z_count+1)=[e_wave3'];
      i_fft_wave3(1:length(t),z_count+1)=[abs(es_wave3').^2];
      i_fftnorm_wave3(1:length(t),z_count+1)=[abs(es_wave3').^2/max(abs(es_wave3').^2)];

      %
      distance(z_count+1)=[zeta];
   end
end;
power_wave1_L(gradual)=(1e+3)*power_wave1(1)/(2*tmax);
power_wave2_L(gradual)=(1e+9)*power_wave2(end)/(2*tmax);
power_wave3_L(gradual)=(1e+3)*power_wave3(1)/(2*tmax);

eta_L(gradual)=100*power_wave2_L(gradual)*(1e-3)/(power_wave3_L(gradual)*power_wave1_L(gradual));
end
%% signal gain
        x = 1e3*L*x;
        gain = 10*log10( power_wave1);
        signal_gain = gain -  10*log10(power_wave1(1));
 %      hold on
        figure(1);
        plot(x, signal_gain);
        xlabel 'Length (mm)'
        ylabel 'Signal Gain (dB)'
        grid on
        grid minor
        fng = sprintf('%0.2f nm, pump wl = %0.2f nm, pump = %0.2f mW Signal gain.eps', wl, pl, pump);
        saveas(gcf, fng);
%% Pump power
%       hold on
        figure(2)
        pump_pwr = power_wave3*f0;
        plot(x, pump_pwr);
        xlabel 'Length (mm)'
        ylabel 'Pump Power (W)'
        grid on
        grid minor
        saveas(gcf, 'pump.eps');
%% signal power
 %      hold on
        figure(3)
        signal_pwr = power_wave1*f0;
        plot(x, signal_pwr);
        xlabel 'Length (mm)'
        ylabel 'Signal Power (W)'
        grid on
        grid minor
        saveas(gcf, 'signal_power.eps');

% %% Total signal and idler gain
%       gain = 10*log10(power_wave2./power_wave2(1) + power_wave1./power_wave1(1)-1);
%       %gain = gain -  10*log10(power_wave1(1));
% %      hold on
%       figure(4)
%       plot(x, gain);
%       xlabel 'Length (mm)'
%       ylabel 'PSA Gain (dB)'
%       grid on
%        grid minor
%        saveas(gcf, 'TotalGain.eps');

%%  Idler gain
        gain = 10*log10( power_wave2);
        idlergain = gain -  10*log10(power_wave2(1));
%     hold on
        figure(5)
        plot(x, idlergain);
        xlabel 'Length (mm)'
        ylabel 'Idler Gain (dB)'
        grid on
        grid minor
        saveas(gcf, 'IdlerGain.eps');
 %%
      rst(:,1) = x;
      rst(:, 2) = signal_gain;
      rst(:, 3) = idlergain;
            rst(:, 4) = pump_pwr;
            rst(:, 5) = signal_pwr;
      fn = sprintf('%0.2f nm, pump wl = %0.2f nm, pump = %0.2f mW PSA gain.txt', wl, pl, pump);
      save (fn, 'rst', '-ascii', '-tabs');

end
