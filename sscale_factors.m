alpha = 0.2;	% power loss (dB/km)
beta2 = 21.67*1e-27;% ps²/km
betap = [0,0, beta2];% dispersion polynomial
gamma = 1.27*1e-3; % Non Linearity coefficient ((W.km)^-1)

L=2000*1e03;% propagation distance

P0=2/(gamma*L); % Power (W)
T0 = sqrt(abs(beta2)*L/2); % Time Period (s or ps ?)
%scaler = 5e9; % original signal has bandwidth -5*10/2:1:5*10/2 where 10 is the normalized spreading and 1/2 is because the time window is 100 (100/2 = 50)

F_carrier = 193.41e12; % Hz
h_planck = 6.626070040e-34; % Js
Kt = 1.13; % phonon occupancy factor at room temperature for Raman amplification in typical optic-fiber
sigma2 = alpha/(10*log10(exp(1))) * h_planck * F_carrier * Kt; % spectral density of noise (eq 56 of the paper)
sigma2n = sigma2 * L / (P0*T0); % normalization of the variance of ASE noise
