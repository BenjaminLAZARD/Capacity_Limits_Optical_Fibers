function u = strangesolitonpulse(t,t0,omega0, A1, A2)

u = A1*sech(A1*(t-t0)).*exp(1i*omega0*t)+ A2*sech(A2*(t+t0)).*exp(-1i*omega0*t);
end

