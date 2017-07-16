function s = WDM_demod (yt,t, W0, N, M)
	% row of s are symbols across channels in frequency
	% columns are symbols in time
	
	
	nt=size(t);
	dt=t(2)-t(1);
	
	p=sinc(W0*t);
	
	%p=conj(p);
	
	Ts=1/W0;

	N_s=floor(Ts/dt);

	u=0;

	n0=(N-1)/2+1;
    m0=(M-1)/2+1;
	
	%demodulation according to the equation (2) of https://arxiv.org/pdf/1302.2875.pdf
	%operated by projections
	for k = -(N-1)/2:(N-1)/2 % frequency
		for l = -(M-1)/2:(M-1)/2 % =over time
			s(n0+k,m0+l) = sum(yt.* circshift(p,[1,l*N_s]).*exp(-1i*2*pi*k*W0*t));
		end
	end
	s=s/norm(p)^2;
end