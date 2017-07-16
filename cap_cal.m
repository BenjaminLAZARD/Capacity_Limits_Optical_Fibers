
[a,b]=size(py_x);

[ rings, phases, P, s,svec ] = constel( Mrings, Nphases, Max_ampl );
cost=abs(svec).^2;

px0=ones(1,a)/a;
e=ones(1,a);

[C,px] = fmincon(@IXY,px0,[],[],e,1,0,1);

Pc2=px*svec';

Pc2_dbm=10*log10(Pc2*P0/1e-3); % dbm
SNR=10*log10(Pc2*P0/( sigma2*Bu*L  ));
