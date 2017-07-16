%check alpha
% check noise power 
%paper: p=0dbm, SNR ~ 20 dB
clc
clear all
close all

h = waitbar(0,'Initializing waitbar...');% create a progress bar
seconds = cputime;

sscale_factors

N_channels = 5; % Number of WDM channels (users)
N_users = 1; %%% Number of symbolsper user

W0 = 1;% per_user bandwidth

%%% TX constellation
Mrings = 2^3;% number of rings (should be a power of 2)
Nphases = 2^3;% number of points on each ring
Max_ampl = 0.01; % in the normalized model
[rings, phases,Pc] = constel( Mrings, Nphases, Max_ampl );

%%% Rx constellation
Mrings_Rx = 2^8;
Nphases_Rx = 2^8;
Max_ampl_Rx = Max_ampl*1.1; % up to 10% noise is the heck of a lot
dring_Rx = Max_ampl_Rx / Mrings_Rx;
dphases_Rx = 2*pi/Nphases_Rx;
[ rings_Rx, phases_Rx,P_Rx] = constel( Mrings_Rx, Nphases_Rx, Max_ampl_Rx );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  time axis
T = 100;% time window (period)
nt = 2^10;% number of points
dt = T/nt;% timestep (dt)
t = -T/2:dt:T/2-dt;% time vector

%%% frequency axis
F=1/dt;%Frequency Window / BandWidth
df=1/T;%Frequency step
f=-F/2:df:F/2-df; %frequencies vector

%%% distance
z=1;% normalized propagation distance
nz=1000; % number of steps
dz=z/nz; % step size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W0 = 1/(floor(1/(dt*W0))*dt);
Bn = N_channels*W0;
Bu=Bn/T0
b2n = sigma2n * Bn * 1/nz;
%b2=sigma2*Bu*L/nz;


Pc_un=10*log10(Pc*P0/1e-3) % dbm
SNR=10*log10(Pc*P0/( sigma2*Bu*L  ))

n0=(N_channels-1)/2+1; % central symbols
m0=(N_users-1)/2+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% starting simulation

N_simu = 10; % then set to 1000
rand('seed'); %

s_hat = zeros(Mrings, N_simu); % estimated symbols for each ring and noise simulation

for ring_index = 1: Mrings % for each ring
	
	s = zeros(N_channels,N_users);% the symbols
    s(n0,m0) = rings(ring_index);% at this point we have 5 chanels modulated over frequency.
    temp= s(n0,m0);
	
	%%%%%%%% For each simulation (to average the results, and get the expectancy)
	for simu_index = 1:N_simu
		
        s = randsymb(rings, phases, N_channels,N_users);
	
        s(n0,m0) = temp;% at this point we have 5 chanels modulated over frequency.		
		
		%%%%%%%% Create the signal 
		u0 = WDM_mod (t, s, W0);% time domain
		
		%% spatial propagation with non linearities (normalized equations)
		%uz = ssprop(transpose(u0),dt,dz,nz, 0, [0 0 -2], 2); % time domain
	    [uz, Pn(ring_index, simu_index)] = sspropn(transpose(u0),dt,dz,nz, 0, [0 0 -2], 2, b2n); % time domain

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

%Calculate noise during propagation
noise = mean(mean(Pn))
SNR_real = 10*log10(Pc/noise)

close(h);% close progress bar

er = plots( u0,uz,uzt_filtered,uzb,t,f);

const_rx_plot( s_hat,rings,phases)

prob = prob_mapping(s_hat, dring_Rx, dphases_Rx);
figure
plot(prob,'.')

py_x = prob_mat(prob, Mrings, Nphases, Mrings_Rx, Nphases_Rx, dring_Rx, dphases_Rx);
plot_proba_matrix(py_x, Mrings, Nphases, Mrings_Rx, Nphases_Rx); %only plots the non circular rotated constellation probabilities

%I = mutual_inf(py_x, Mrings, Nphases, Mrings_Rx, Nphases_Rx)

px = ones(1, Mrings*Nphases);
for i=1:Nphases*Mrings
    px(i) = 1/(Mrings*Nphases);
end

I = IXY(py_x, px)

