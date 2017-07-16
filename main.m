clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% defining time, frequency and distance parameters
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

%%% Number of users to filter (in case of frequency multiplexing)
N_users = 1;

%%% distance
z=1;% normalized propagation distance
nz=1000; % number of steps
dz=z/nz; % step size

%%%cf sscale_factors.m to get back to orignal units
scaler = 1;
tt = t/scaler;
ff = f*scaler;
alpha = 0.1;

%B=0.5;
%s=1;
%u0 = s*sinc(B*t);%raisedCosinus(t, 0.2, 1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% defining the input
%%% modulation symbols as per time and frequency
s=0.3*[0.5 1 0.3+0.4j 0.4 0.9]';
[N,M]=size(s);

%%% modulation
u0 = WDM_mod (t, s, W0);% time domain
u0f=fftshift(fft(u0));% frequency domain

%%% spatial propagation with non linearities (normalized equations)
%uz = sspropn(transpose(u0),dt,dz,nz, alpha, [0 0 -2], 2, 100); % time domain with noise
uz = ssprop(transpose(u0),dt,dz,nz, alpha, [0 0 -2], 2); % time domain
uzf = fftshift(fft(uz)); % frequency domain

%%% filter the Canal Of Interest (specific user for example)
[u0t_filtered, u0f_filtered]=Filter(u0, df, N_users*W0);
[uzt_filtered, uzf_filtered]=Filter(uz, df, N_users*W0);

%%% back_propagation (what was the signal at the beginning, knowing the parameters of the fiber ?)
uzb = ssprop(uzt_filtered,dt,dz,nz, -alpha, [0 0 2],-2);
uzb = transpose(uzb); % time domain
uzbf = fftshift(fft(uzb)); % frequency domain

%%% estimator of the original symbols by demodulation
shat = WDM_demod (uzb, t, W0, N, M); % symbols estimated




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% comparing input and output
%%% plot in time domain
figure(1)
plot(tt, abs(u0), "-b;modulated signal;",
	 tt, abs(uz),"-r;spatial propagation;",
	 tt, abs(uzt_filtered),"-g;selecting user = filtering;",
	 tt, abs(uzb), "-k;backpropagation;"
	)
xlabel ("time (s)");
ylabel ("Power (W)");
title ("Comparing the signal after different steps in time domain");

%%% plot in frequency domain
figure(2)
plot(ff, abs(u0f), "-b;modulated signal;",
	 ff, abs(uzf), "-r;spatial propagation;",
	 ff, abs(uzf_filtered), "-g;selecting user = filtering;",
	 ff, abs(uzbf), "-k;backpropagation;"
	)
xlabel ("frequency (Hz)");
ylabel ("Power (W)");
title ("Comparing the signal after different steps in frequency domain");

%%% evaluate the error
%et = norm(uzb-u0t_filtered)/norm(u0t_filtered) % very peaky
ef = norm(uzbf-u0f_filtered)/norm(u0f_filtered)