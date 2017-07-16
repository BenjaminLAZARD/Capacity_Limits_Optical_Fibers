function const_rx_plot( shat, rings,phases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% comparing input and output
%%% plot in time domain

rmax=rings(end);

[N,M]=size(shat);

Colors = {'r','g','b','m','c','k','r','g','b','m','c','k','r','g','b','m','c','k'};

figure

for k=1:N
    
plot(real(shat(k,:)), imag(shat(k,:)),'.','Color',Colors{k})
xlabel ('Re(s)');
ylabel ('Im(s)');
title ('constellation, RX');
hold on

th = 0:pi/50:2*pi;
xunit = rings(k) * cos(th);
yunit = rings(k) * sin(th);
plot(xunit, yunit,'Color',Colors{k});
hold on

axis(1.2*[-rmax rmax -rmax rmax])
end


%function h = circle(x,y,r)
%hold on
% for i=1:length(rings)
%     th = 0:pi/50:2*pi;
%     xunit = rings(i) * cos(th);
%     yunit = rings(i) * sin(th);
%     plot(xunit, yunit);
%     hold on
% end

hold off

end