function [e_wave1,e_wave2,e_wave3,e_wave4,es_wave1,es_wave2,es_wave3,es_wave4] = step_4wave_4wave(e_wave1,e_wave2,e_wave3,e_wave4,deltat,deltaz,deltak,gvd_wave1,gvd_wave2,gvd_wave3,gvd_wave4,c_chi3_wave1,c_chi3_wave2,c_chi3_wave3,c_chi3_wave4,gvm_wave1,gvm_wave2,gvm_wave3,gvm_wave4,pphase,alpha_wave1,alpha_wave2,alpha_wave3,alpha_wave4,c_n2a2_wave1,c_n2a2_wave2,c_n2a2_wave3,c_n2a2_wave4)
%% Wave assignment
% 3 and 4 - pump
% 1 - signal
% 2 - idler
%% 1st Linear Half Step
%  wave1: FFT and linear half step
    j = sqrt(-1);
    es_wave1=deltat*fftshift(fft(e_wave1)); % FFT
    fact_wave1=-j.*gvd_wave1.*deltaz-j.*gvm_wave1.*deltaz-0.5*alpha_wave1.*deltaz;
    es_wave1=exp(fact_wave1*0.5).*es_wave1;
%  FF: IFFT
    e_wave1=(1/deltat)*ifft(ifftshift(es_wave1));
    
%  wave2: FFT and linear half step
    es_wave2=deltat*fftshift(fft(e_wave2)); % FFT
    fact_wave2=-j.*gvd_wave2.*deltaz-j.*gvm_wave2.*deltaz-0.5*alpha_wave2.*deltaz;
    es_wave2=exp(fact_wave2*0.5).*es_wave2;
%  FF: IFFT
    e_wave2=(1/deltat)*ifft(ifftshift(es_wave2));
    
%  wave3: FFT and linear half step
    es_wave3=deltat*fftshift(fft(e_wave3)); % FFT
    fact_wave3=-j.*gvd_wave3.*deltaz-j.*gvm_wave3.*deltaz-0.5*alpha_wave3.*deltaz;
    es_wave3=exp(fact_wave3*0.5).*es_wave3;
% SH: IFFT
    e_wave3=(1/deltat)*ifft(ifftshift(es_wave3));
    
%   wave4: FFT and linear half step
    es_wave4=deltat*fftshift(fft(e_wave4)); % FFT
    fact_wave4=-j.*gvd_wave4.*deltaz-j.*gvm_wave4.*deltaz-0.5*alpha_wave4.*deltaz;
    es_wave4=exp(fact_wave4*0.5).*es_wave4;
%  FF: IFFT
    e_wave4=(1/deltat)*ifft(ifftshift(es_wave4));

%%  quadratic nonlinear step
% Pump 1: m 
% Pump 2: l
% k1 -> A2
% j1 -> A3 
    % Pump 1
    m1 = c_chi3_wave3*conj(e_wave4).*(e_wave2.*e_wave1).*exp(-j*pphase) + c_n2a2_wave3*(conj(e_wave3).*e_wave3).*e_wave3;
    % Pump 2
    l1 = c_chi3_wave4*conj(e_wave3).*(e_wave2.*e_wave1).*exp(-j*pphase) + c_n2a2_wave4*(conj(e_wave4).*e_wave4).*e_wave4;
    % Idler
    k1 = c_chi3_wave2*(e_wave4.*e_wave3).*conj(e_wave1).*exp(-j*pphase) + c_n2a2_wave2*(conj(e_wave2).*e_wave2).*e_wave2;
    % Signal
    j1 = c_chi3_wave1*(e_wave4.*e_wave3).*conj(e_wave2).*exp(-j*pphase) + c_n2a2_wave1*(conj(e_wave1).*e_wave1).*e_wave1;
    
    % Pump 1
    m2 = c_chi3_wave3*conj(e_wave4+deltaz*l1).*(e_wave2+deltaz*k1).*(e_wave1+deltaz*j1).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave3*conj(e_wave3+deltaz*m1).*(e_wave3+deltaz*m1).*(e_wave3+deltaz*m1);
    % Pump 2
    l2 = c_chi3_wave4*conj(e_wave3+deltaz*m1).*(e_wave2+deltaz*k1).*(e_wave1+deltaz*j1).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave4*conj(e_wave4+deltaz*l1).*(e_wave4+deltaz*l1).*(e_wave4+deltaz*l1);
    % Idler
    k2 = c_chi3_wave2*(e_wave4+deltaz*l1).*(e_wave3+deltaz*m1).*conj(e_wave1+deltaz*j1).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave2*conj(e_wave2+deltaz*k1).*(e_wave2+deltaz*k1).*(e_wave2+deltaz*k1);
    % Signal
    j2 = c_chi3_wave1*(e_wave4+deltaz*l1).*(e_wave3+deltaz*m1).*conj(e_wave2+deltaz*k1).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave1*conj(e_wave1+deltaz*j1).*(e_wave1+deltaz*j1).*(e_wave1+deltaz*j1);
    
    % Pump 1
    m3 = c_chi3_wave3*conj(e_wave4+deltaz*l2).*(e_wave2+deltaz*k2).*(e_wave1+deltaz*j2).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave3*conj(e_wave3+deltaz*m2).*(e_wave3+deltaz*m2).*(e_wave3+deltaz*m2);
    % Pump 2
    l3 = c_chi3_wave4*conj(e_wave3+deltaz*m2).*(e_wave2+deltaz*k2).*(e_wave1+deltaz*j2).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave4*conj(e_wave4+deltaz*l2).*(e_wave4+deltaz*l2).*(e_wave4+deltaz*l2);
    % Idler
    k3 = c_chi3_wave2*(e_wave4+deltaz*l2).*(e_wave3+deltaz*m2).*conj(e_wave1+deltaz*j2).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave2*conj(e_wave2+deltaz*k2).*(e_wave2+deltaz*k2).*(e_wave2+deltaz*k2);
    % Signal
    j3 = c_chi3_wave1*(e_wave4+deltaz*l2).*(e_wave3+deltaz*m2).*conj(e_wave2+deltaz*k2).*exp(-j*(pphase+0.5*deltak*deltaz)) + c_n2a2_wave1*conj(e_wave1+deltaz*j2).*(e_wave1+deltaz*j2).*(e_wave1+deltaz*j2);
    
    % Pump 1
    m4 = c_chi3_wave3*conj(e_wave4+deltaz*l3).*(e_wave2+deltaz*k3).*(e_wave1+deltaz*j3).*exp(-j*(pphase+deltak*deltaz)) + c_n2a2_wave3*conj(e_wave3+deltaz*m3).*(e_wave3+deltaz*m3).*(e_wave3+deltaz*m3);
    % Pump 2
    l4 = c_chi3_wave4*conj(e_wave3+deltaz*m3).*(e_wave2+deltaz*k3).*(e_wave1+deltaz*j3).*exp(-j*(pphase+deltak*deltaz)) + c_n2a2_wave4*conj(e_wave4+deltaz*l3).*(e_wave4+deltaz*l3).*(e_wave4+deltaz*l3);
    % Idler
    k4 = c_chi3_wave2*(e_wave4+deltaz*l3).*(e_wave3+deltaz*m3).*conj(e_wave1+deltaz*j3).*exp(-j*(pphase+deltak*deltaz)) + c_n2a2_wave2*conj(e_wave2+deltaz*k3).*(e_wave2+deltaz*k3).*(e_wave2+deltaz*k3);
    % Signal
    j4 = c_chi3_wave1*(e_wave4+deltaz*l3).*(e_wave3+deltaz*m3).*conj(e_wave2+deltaz*k3).*exp(-j*(pphase+deltak*deltaz)) + c_n2a2_wave1*conj(e_wave1+deltaz*j3).*(e_wave1+deltaz*j3).*(e_wave1+deltaz*j3);
    
    e_wave1=e_wave1+(j1+2*j2+2*j3+j4)*deltaz/6;
    e_wave2=e_wave2+(k1+2*k2+2*k3+k4)*deltaz/6;
    e_wave3=e_wave3+(m1+2*m2+2*m3+m4)*deltaz/6;
    e_wave4=e_wave4+(l1+2*l2+2*l3+l4)*deltaz/6;

%% 2nd Linear Half Step   
%  wave1: FFT and linear half step
    es_wave1=deltat*fftshift(fft(e_wave1)); % FFT
    fact_wave1=-j.*gvd_wave1.*deltaz-j.*gvm_wave1.*deltaz-0.5*alpha_wave1.*deltaz;
    es_wave1=exp(fact_wave1*0.5).*es_wave1;
%  FF: IFFT
    e_wave1=(1/deltat)*ifft(ifftshift(es_wave1));
    
%  wave2: FFT and linear half tep
    es_wave2=deltat*fftshift(fft(e_wave2)); % FFT
    fact_wave2=-j.*gvd_wave2.*deltaz-j.*gvm_wave2.*deltaz-0.5*alpha_wave2.*deltaz;
    es_wave2=exp(fact_wave2*0.5).*es_wave2;
%  FF: IFFT
    e_wave2=(1/deltat)*ifft(ifftshift(es_wave2));
    
%  wave3: FFT and linear half step
    es_wave3=deltat*fftshift(fft(e_wave3)); % FFT
    fact_wave3=-j.*gvd_wave3.*deltaz-j.*gvm_wave3.*deltaz-0.5*alpha_wave3.*deltaz;
    es_wave3=exp(fact_wave3*0.5).*es_wave3;
% SH: IFFT
    e_wave3=(1/deltat)*ifft(ifftshift(es_wave3));

%  wave4: FFT and linear half step
    es_wave4=deltat*fftshift(fft(e_wave4)); % FFT
    fact_wave4=-j.*gvd_wave4.*deltaz-j.*gvm_wave4.*deltaz-0.5*alpha_wave4.*deltaz;
    es_wave4=exp(fact_wave4*0.5).*es_wave4;
% SH: IFFT
    e_wave4=(1/deltat)*ifft(ifftshift(es_wave4));
