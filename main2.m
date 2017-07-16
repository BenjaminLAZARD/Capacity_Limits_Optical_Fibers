clc
clear all
close all

h = waitbar(0,'Initializing waitbar...');% create a progress bar
seconds = cputime;

sscale_factors

N_channels = 5; % Number of WDM channels (users)

N_users = 1; %%% Number of symbolsper user


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% defining time, frequency, distance, normalization, and simulation parameters
%%% time
T = 100;% time window (period)
nt = 2^10;% number of points
dt = T/nt;% timestep (dt)
t = -T/2:dt:T/2-dt;% time vector

%%% frequency
F=1/dt;%Frequency Window / BandWidth
df=1/T;%Frequency step
f=-F/2:df:F/2-df; %frequencies vector
W0 = 1;%normalized frequency of the modulating functions
W0 = 1/(floor(1/(dt*W0))*dt);

%%% distance
z=1;% normalized propagation distance
nz=1000; % number of steps
dz=z/nz; % step size

%%% Number of simulations
N_simu = 10; % then set to 1000
rand('seed'); % setting a seed (is it the right way to do so ?)

% 
 n0=(N_channels-1)/2+1; % central symbols
 m0=(N_users-1)/2+1;
% %%% normalization cf sscale_factors.m to get back to orignal units
% %scaler = 5e9;
% alpha = 0.2;	% power loss (dB/km)
% beta2 = 21.67*1e-27;% ps²/km
% betap = [0,0, beta2];% dispersion polynomial
% gamma = 1.27*1e-3; % Non Linearity coefficient ((W.km)^-1)
% L = 2000*1e03;% propagation distance
% P0 = 2/(gamma*L); % Power (W)
% T0 = sqrt(abs(beta2)*L/2); % Time Period (s)
% F_carrier = 193.41e12; % Hz
% h_planck = 6.626070040e-34; % Js
% Kt = 1.13; % phonon occupancy factor at room temperature for Raman amplification in typical optic-fiber
% sigma2 = alpha/(10*log10(exp(1)))* h_planck * F_carrier * Kt; % spectral density of noise (eq 56 of the paper)
% sigma2n = sigma2 * L / (P0*T0); % normalization of the variance of ASE noise
% b2 = sigma2n * Bn * 1/nz;

Bn = N_channels*W0;



%%% Ring constellations parameters
Mrings = 2^3;% number of rings (should be a power of 2)
Nphases = Mrings;% number of points on each ring
Max_ampl = 1.5; % in the normalized model
% dring = Max_ampl/Mrings; %  uniformly distributed rings
% dphases = 2*pi/Nphases; % uniformlfy distributed phases
% %%% Generate the Tx Ring constellation
% rings = dring:dring:Max_ampl;
% phases = 0:dphases:2*pi-dphases; % uniformly distributed phases
% constellation = zeros(Mrings, Nphases);
% for i = 1:Mrings
% 	constellation(i,:) = rings(i)*exp(1i*phases);
% 	figure(1)
% 	plot(constellation(i,:), 'o');
% 	hold on;
% end

 [rings, phases] = constel( Mrings, Nphases, Max_ampl );


%%% Rx constellation parameters (a finner and slightly larger mapping of |R²)
Mrings_Rx = 2^8;
Nphases_Rx = 2^8;
Max_ampl_Rx = Max_ampl*1.1; % up to 10% noise is the heck of a lot
% dring_Rx = Max_ampl_Rx / Mrings_Rx; % uniformly distributed rings
% dphases_Rx = 2*pi/Nphases_Rx; % uniformlfy distributed phases
% %%% Generate the Rx Ring constellation
% rings_Rx = dring_Rx:dring_Rx:Max_ampl_Rx;
% phases_Rx = 0:dphases_Rx:2*pi-dphases_Rx; % uniformly distributed phases
% constellation_Rx = zeros(Mrings_Rx, Nphases_Rx);
% for i = 1:Mrings_Rx
% 	constellation_Rx(i,:) = rings_Rx(i)*exp(1i*phases_Rx);
% %	figure(2)
% %	plot(constellation_Rx(i,:), 'o');
% %	hold on;
% end

[ rings_Rx, phases_Rx,Pc] = constel( Mrings_Rx, Nphases_Rx, Max_ampl_Rx );

Pc_un=10*log10(Pc*P0/1e-3);

SNR=10*log10(Pc/( b2*nz  ))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% starting simulation

s_hat = zeros(Mrings, N_simu); % estimated symbols for each ring and noise simulation

for ring_index = 4: 4% Mrings % for each ring
	
	s = zeros(N_channels,N_users);% the symbols
	%[N, M] = size(s);
   
	
	%%%%%%%% For each simulation (to average the results, and get the expectancy)
	for simu_index = 1:N_simu
		
		%%%%%%%% Select the symbols to be sent.
		%% generate noisy values for the symbols from the constellation:
% 		si = round(unifrnd(1, Mrings, N_channels, 2));%pick 5 indexes in [1,Mrings]²
% 		for i = 1:N_channels %assign the corresponding constellation values to s
% 			s(:,i) = constellation(si(i,1), si(i,2));
% 		end

        s = randsymb(rings, phases, N_channels,N_users);
		%% current user = current ring: set the middle symbol to the current ring radius
		
         s(n0,m0) = rings(ring_index);% at this point we have 5 chanels modulated over frequency.
         s0=s;
		
		
		%%%%%%%% Create the signal
        
		u0 = WDM_mod (t, s, W0);% time domain
		
		%% spatial propagation with non linearities (normalized equations)
		%uz = ssprop(transpose(u0),dt,dz,nz, 0, [0 0 -2], 2); % time domain
	    uz = sspropn(transpose(u0),dt,dz,nz, 0, [0 0 -2], 2,b2); % time domain

		%% filter the Canal Of Interest (specific user for example)
		[uzt_filtered, uzf_filtered] = Filter(uz, df, W0); % this is equivalent ot ROADM (cf.paper)
        %uzt_filtered = uz;
        %uzf_filtered = fft(fftshift(uz));
		
		%% back_propagation (what was the signal at the beginning, knowing the parameters of the fiber ?)
		uzb = transpose(ssprop(uzt_filtered,dt,dz,nz, 0, [0 0 2], -2)); % time domain
		
		%% demodulation
        s_all = WDM_demod(uzb, t, W0, N_channels, N_users);
		s_hat(ring_index, simu_index) = s_all(n0,m0); % estimated central symbols
		
		%% printing rate of simulation
		%printf("%d out of %d \n", ring_index*simu_index, Mrings * N_simu);%printf("ring %d/%d, simu %d/%d \n", ring_index, Mrings, simu_index, N_simu);
		perc = (N_simu*(ring_index-1) + simu_index) / (Mrings * N_simu);
		waitbar(perc, h, sprintf('%d%%', perc*100))
		sprintf('time elapsed:%s s', cputime - seconds);
	end
end

close(h);% close progress bar


figure
const_rx_plot( s_hat,rings,phases)

er = plots( u0,uz,uzt_filtered,uzb,t,f);

%%%%%%%% Output probabilities for each value of Chanel of Interest symbol.
%prob = prob_mapping(s_hat, dring_Rx, dphases_Rx);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plot
% figure(3)
% for i = 1:Mrings
% 	plot(prob(i,:), 'o')
%     hold on
% end
%close all

