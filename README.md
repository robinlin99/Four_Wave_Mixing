# Four_Wave_Mixing

**Four_Wave_Mixing** is a numerical simulation for Chi-3 degenerate and nondegenerate four-wave mixing (FWM) process. Currently, the nondegenerate FWM solver is still being developed. 

Two main files will be used in the simulation:
1. **PIA_Gain.m** is the main program. The function takes inputs **wl**, **pump**, and **pl**. These correspond to signal wavelength, pump power, and pump wavelength, respectively. As output, idler gain and power, the signal gain and power, and pump power, as functions of waveguide propagation length, will be generated. 
2. **step_4wave.m** is the numerical solver for FWM. It takes into account fiber disperion, two-photon absorption, and self-phase modulation. Split-step fourier method is used for the linear and nonlinear parts of the differential equation. The quadratic nonlinear step utilizes a fourth-order runge-kutta algorithm. 
