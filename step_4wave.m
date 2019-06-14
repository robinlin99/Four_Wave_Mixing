function [e_wave1,e_wave2,e_wave3,es_wave1,es_wave2,es_wave3] = step_4wave(e_wave1,e_wave2,e_wave3,deltat,deltaz,deltak,gvd_wave1,gvd_wave2,gvd_wave3,c_chi3_wave1,c_chi3_wave2,c_chi3_wave3,gvm_wave1,gvm_wave2,gvm_wave3,pphase,alpha_wave1,alpha_wave2,alpha_wave3,c_n2a2_wave1,c_n2a2_wave2,c_n2a2_wave3)

%% 1st Linear Half Step
%  wave1: FFT and linear half step
    j = sqrt(-1);
    es_wave1=deltat*fftshift(fft(e_wave1)); % FFT
    fact_wave1=-j.*gvd_wave1.*deltaz+j.*gvm_wave1.*deltaz-0.5*alpha_wave1.*deltaz;
    es_wave1=exp(fact_wave1*0.5).*es_wave1;
%  FF: IFFT
    e_wave1=(1/deltat)*ifft(ifftshift(es_wave1));
    
%  wave2: FFT and linear half tep4
    es_wave2=deltat*fftshift(fft(e_wave2)); % FFT
    fact_wave2=-j.*gvd_wave2.*deltaz+j.*gvm_wave2.*deltaz-0.5*alpha_wave2.*deltaz;
    es_wave2=exp(fact_wave2*0.5).*es_wave2;
%  FF: IFFT
    e_wave2=(1/deltat)*ifft(ifftshift(es_wave2));
    
%  wave3: FFT and linear half step
    es_wave3=deltat*fftshift(fft(e_wave3)); % FFT
    fact_wave3=-j.*gvd_wave3.*deltaz+j.*gvm_wave3.*deltaz-0.5*alpha_wave3.*deltaz;
    es_wave3=exp(fact_wave3*0.5).*es_wave3;
% SH: IFFT
    e_wave3=(1/deltat)*ifft(ifftshift(es_wave3));


%%  quadratic nonlinear step
%  two-variable fourth-order Runge-Kutta algorithm
% m1 -> pump (degenerate case) -> e_wave1
% k1 -> A2 -> e_wave2
% j1 -> A3 -> e_wave3
    % A3
    j1=c_chi3_wave1*(e_wave3.*e_wave3).*conj(e_wave2).*exp(-j*pphase) + c_n2a2_wave1*(conj(e_wave1).*e_wave1).*e_wave1;
    % A2
    k1=c_chi3_wave2*(e_wave3.*e_wave3).*conj(e_wave1).*exp(-j*pphase) + c_n2a2_wave2*(conj(e_wave2).*e_wave2).*e_wave2;
    % Pump
    m1=c_chi3_wave3*conj(e_wave3).*(e_wave2.*e_wave1).*exp(j*pphase) + c_n2a2_wave3*(conj(e_wave3).*e_wave3).*e_wave3;
    

    j2=c_chi3_wave1*(e_wave3 + deltaz*m1).*(e_wave3 + deltaz*m1).*conj(e_wave2 + deltaz*k1).*exp(-j*(pphase+0.5*deltak*deltaz))+ c_n2a2_wave1*(conj(e_wave1+0.5*deltaz*j1).*(e_wave1+0.5*deltaz*j1)).*(e_wave1+0.5*deltaz*j1);
    k2=c_chi3_wave2*(e_wave3 + deltaz*m1).*(e_wave3 + deltaz*m1).*conj(e_wave1 + deltaz*j1).*exp(-j*(pphase+0.5*deltak*deltaz))+ c_n2a2_wave2*(conj(e_wave2+0.5*deltaz*k1).*(e_wave2+0.5*deltaz*k1)).*(e_wave2+0.5*deltaz*k1);
    m2=c_chi3_wave3*conj(e_wave3 + deltaz*m1).*(e_wave2 + deltaz*k1).*(e_wave1 + deltaz*j1).*exp(j*(pphase+0.5*deltak*deltaz))+ c_n2a2_wave3*(conj(e_wave3+0.5*deltaz*m1).*(e_wave3+0.5*deltaz*m1)).*(e_wave3+0.5*deltaz*m1);


    j3=c_chi3_wave1*(e_wave3 + deltaz*m2).*(e_wave3 + deltaz*m2).*conj(e_wave2 + deltaz*k2).*exp(-j*(pphase+0.5*deltak*deltaz))+ c_n2a2_wave1*(conj(e_wave1+0.5*deltaz*j2).*(e_wave1+0.5*deltaz*j2)).*(e_wave1+0.5*deltaz*j2);
    k3=c_chi3_wave2*(e_wave3 + deltaz*m2).*(e_wave3 + deltaz*m2).*conj(e_wave1 + deltaz*j2).*exp(-j*(pphase+0.5*deltak*deltaz))+ c_n2a2_wave2*(conj(e_wave2+0.5*deltaz*k2).*(e_wave2+0.5*deltaz*k2)).*(e_wave2+0.5*deltaz*k2);
    m3=c_chi3_wave3*conj(e_wave3 + deltaz*m2).*(e_wave2 + deltaz*k2).*(e_wave1 + deltaz*j2).*exp(j*(pphase+0.5*deltak*deltaz))+ c_n2a2_wave3*(conj(e_wave3+0.5*deltaz*m2).*(e_wave3+0.5*deltaz*m2)).*(e_wave3+0.5*deltaz*m2);

    
    j4=c_chi3_wave1*(e_wave3 + deltaz*m3).*(e_wave3 + deltaz*m3).*conj(e_wave2 + deltaz*k3).*exp(-j*(pphase+deltak*deltaz))+ c_n2a2_wave1*(conj(e_wave1+deltaz*j3).*(e_wave1+deltaz*j3)).*(e_wave1+deltaz*j3);
    k4=c_chi3_wave2*(e_wave3 + deltaz*m3).*(e_wave3 + deltaz*m3).*conj(e_wave1 + deltaz*j3).*exp(-j*(pphase+deltak*deltaz))+ c_n2a2_wave2*(conj(e_wave2+deltaz*k3).*(e_wave2+deltaz*k3)).*(e_wave2+deltaz*k3);
    m4=c_chi3_wave3*conj(e_wave3 + deltaz*m3).*(e_wave2 + deltaz*k3).*(e_wave1 + deltaz*j3).*exp(j*(pphase+deltak*deltaz))+ c_n2a2_wave3*(conj(e_wave3+deltaz*m3).*(e_wave3+deltaz*m3)).*(e_wave3+deltaz*m3);
    

    e_wave1=e_wave1+(j1+2*j2+2*j3+j4)*deltaz/6;
    e_wave2=e_wave2+(k1+2*k2+2*k3+k4)*deltaz/6;
    e_wave3=e_wave3+(m1+2*m2+2*m3+m4)*deltaz/6;


%% 2nd Linear Half Step   
%  wave1: FFT and linear half step
    es_wave1=deltat*fftshift(fft(e_wave1)); % FFT
    fact_wave1=-j.*gvd_wave1.*deltaz+j.*gvm_wave1.*deltaz-0.5*alpha_wave1.*deltaz;
    es_wave1=exp(fact_wave1*0.5).*es_wave1;
%  FF: IFFT
    e_wave1=(1/deltat)*ifft(ifftshift(es_wave1));
    
%  wave2: FFT and linear half tep
    es_wave2=deltat*fftshift(fft(e_wave2)); % FFT
    fact_wave2=-j.*gvd_wave2.*deltaz+j.*gvm_wave2.*deltaz-0.5*alpha_wave2.*deltaz;
    es_wave2=exp(fact_wave2*0.5).*es_wave2;
%  FF: IFFT
    e_wave2=(1/deltat)*ifft(ifftshift(es_wave2));
    
%  wave3: FFT and linear half step
    es_wave3=deltat*fftshift(fft(e_wave3)); % FFT
    fact_wave3=-j.*gvd_wave3.*deltaz+j.*gvm_wave3.*deltaz-0.5*alpha_wave3.*deltaz;
    es_wave3=exp(fact_wave3*0.5).*es_wave3;
% SH: IFFT
    e_wave3=(1/deltat)*ifft(ifftshift(es_wave3));










