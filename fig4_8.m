T = 100;                                % time window (period)
nt = 2^12;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector
w = wspace(T,nt);                       % angular frequency vector
vs = fftshift(w/(2*pi));                % shifted for plotting

z = 5;                             		% total distance
nz = 50;                                % total number of steps
nplot = 5;								% number of plots to make
n1 = round(nz/nplot);					% number of steps per plot
nz = n1*nplot;							% total number of steps (revised)
dz = z/nz;                              % step-size

betap = [0,0,-1];						% dispersion polynomial

u0 = gaussian(t,0,2*sqrt(log(2)));		% u0 = exp(-t.^2/2);

zv = (z/nplot)*(0:nplot);
u = zeros(length(t),length(zv));
U = zeros(length(t),length(zv));

u(:,1) = u0;
U(:,1) = fftshift(abs(dt*fft(u(:,1))/sqrt(2*pi)).^2);

for ii = 1:nplot,
  u(:,ii+1) = ssprop(u(:,ii),dt,dz,n1,0,betap,1);
  U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
end

subplot(211);
mesh(zv,t,abs(u).^2,...
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
