f = 193.41*10^12;
T = 1/f;                                % time window (period)
nt = 2^9;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = 500;                      		% total distance
nz = 50000;                              % total number of steps
nplot = 10;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

alpha = 0.2;
beta2 = 21.67;                          
betap = [0,0,beta2];				% dispersion polynomial
gamma = 1.27;

u0 = sinc(2*pi*f*t);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = u0;
U(:,1) = fftshift(abs(dt*fft(u(:,1))/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = ssprop(u(:,ii),dt,dz,n1,alpha,betap,gamma);
  U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
end

figure();

subplot(211);
mesh(zv/(pi/2),t,abs(u).^2, 'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir', 'reverse');
%set(gca,'ZScale',{'log'});
hidden off;
xlim([0 300]);
ylim([-T T]);
xlabel ('Z/Z_0');
ylabel ('(t-\beta_1z)/T_0');
zlabel ('|u(z,t)|^2/P_0');

subplot(212);
mesh(zv/(pi/2),vs,U, ...
     'MeshStyle', 'col', 'EdgeColor', 'black');
set(gca,'YDir', 'reverse');
hidden off;
xlim([0 3]);
ylim([-f f]);
xlabel ('Z/L_d');
ylabel ('(\nu-\nu_0) T_0');
zlabel ('|U(z,\nu)|^2/P_0');