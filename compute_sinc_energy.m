T = 100;% time window (period)
nt = 2^10;% number of points
dt = T/nt;% timestep (dt)
t = -T/2:dt:T/2-dt;% time vector

for i = 1:10:500
	mysinc = sinc(i*t);
	power = mean(abs(mysinc).^2);
	plot(i, power,'bo');
	hold on;
end