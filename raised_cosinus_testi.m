f = 193.41*(10^12) 						%frequency in THz
T = 100;                                % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = 5;                             		% total distance (km)
nz = 1000;                              % total number of steps (1 every m)
nplot = 5;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

alpha = 0.2								% power loss (dB/km)
beta2 = 21.67;							% ps²/km
betap = [0,0, beta2];					% dispersion polynomial
gamma = 1.27							% Non Linearity coefficient ((W.km)^-1)

rho = 0.1									%raisedCosinus parameter
u0 = raisedCosinus(t, rho, 1);		

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = u0;
U(:,1) = fftshift(abs(dt*fft(u(:,1))/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = ssprop(u(:,ii),dt,dz,n1, alpha,betap,gamma);
  U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
end

figure()
subplot(211);
mesh(zv,t,abs(u).^2, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir','reverse');
hidden off;
xlim([0 5]);
ylim([-12 12]);
xlabel ('Z/L_d');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');

subplot(212);
mesh(zv,vs,U,...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir','reverse');
hidden off;
xlim([0 5]);
ylim([-3 3]);
xlabel ('Z/L_d');
ylabel ('(\nu-\nu_0) T_0');
zlabel ('|U(z,\nu)|^2/P_0');