function [ ef ] = plots( u0t,uzt,uzfilt,uzbt,tt,ff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% comparing input and output
%%% plot in time domain

u0f = fftshift(fft(u0t)); % frequency domain
uzf = fftshift(fft(uzt)); % frequency domain
ufilf = fftshift(fft(uzfilt)); % frequency domain
uzbf = fftshift(fft(uzbt)); % frequency domain



figure(1)
plot(tt, abs(u0t), '-b', tt, abs(uzt),'-r', tt, abs(uzfilt),'-g',tt, abs(uzbt), '-k')
xlabel ('time (s)');
ylabel ('|u(t)|');
title ('Signals in time');
legend('input','output','output+filter','output+filter+BP')

%%% plot in frequency domain
figure(2)
plot(ff, abs(u0f), '-b',ff, abs(uzf), '-r',ff, abs(ufilf), '-g',ff, abs(uzbf), '-k')
xlabel ('frequency (Hz)');
ylabel ('|u(f)|');
title ('Signals in the frequency');
legend('input','output','output+filter','output+filter+BP')

%%% evaluate the error
%et = norm(uzb-u0t_filtered)/norm(u0t_filtered) % very peaky
ef = norm(uzbf-u0f)/norm(u0f)

end

