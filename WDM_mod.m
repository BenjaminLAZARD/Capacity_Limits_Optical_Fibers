function u = WDM_mod (t, s, W0)
	
	% s are the symbols sent
	% row of s are symbols across channels in frequency
	% columns are symbols in time
	% W0 is the frequency of the carrier, always a sinc here
	
	[N M] = size(s);
	nt=size(t);
	dt=t(2)-t(1);
	
	p=sinc(W0*t); %carrier
	Ts=1/W0; %Bandwidth
	
	N_s=floor(Ts/dt);%phase shift step

	u=0;

	n0=(N-1)/2+1;
    m0=(M-1)/2+1;
	
	%modulation according to the equation (2) of https://arxiv.org/pdf/1302.2875.pdf
	for k = -(N-1)/2:(N-1)/2 % frequency
		for l =  -(M-1)/2:(M-1)/2 % =over time
			u=u+ s(n0+k,l+m0) * circshift(p,[1 l*N_s]).*exp(1i*2*pi*k*W0*t);
		end
	end
end