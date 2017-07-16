clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

p=sinc(W0*t);



s=[ 1+j 1-2.5j 0.7j; 1.1+j 2.2-j 1.6-1.1*j  ];
[N,M]=size(s);

%s=1.7-j;
% 
% % no propagation
% u = WDM_mod (t, s, W0);
% [N,M] = size(s);
% shat = WDM_demod (u, t, W0, N, M);

% 
% with propagation

z=1;
nz=1000;
dz=z/nz;


u0 = WDM_mod (t, s, W0);
uz = ssprop(transpose(u0),dt,dz,nz, 0, [0 0 -2], 2);
uzb = transpose(ssprop(uz,dt,dz,nz, 0, [0 0 2], -2));
shat = WDM_demod (uzb, t, W0, N, M);



% %%%% check directly
% Ts=1/W0; %Bandwidth
% N_s=floor(Ts/dt);%phase shift step
% pm1 = circshift(p,[1,-1*N_s]);
% p0 = circshift(p,[1,0*N_s]);
% pp1 = circshift(p,[1,1*N_s]);
% 
% sm1=sum(u.*pm1)/norm(p)^2;
% s0=sum(u.*p)/norm(p)^2;
% sp1=sum(u.*pp1)/norm(p)^2;


%%%%%%%%%%%

%uf=fftshift(fft(u));

% figure
% plot(f, abs(uf),'r')
% figure
% plot(t,abs(u),'b');


