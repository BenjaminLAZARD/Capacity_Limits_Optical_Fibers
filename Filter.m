function  [yt,yf]= Filter(ut,df,B)
% filters the time signal ut in bandwidth B with frequency step df
% returns the filtered signal in both time and frequency domains
	
N=numel(ut);
uf=fftshift(fft(ut));

Bn=floor(0.5*B/df);
N0=N/2+1;
I=N0-Bn:N0+Bn;

yf=zeros(1,N);
yf(I)=uf(I);

%yt=ifftshift(ifft(uf));
yt = transpose(ifft(fftshift(yf)));
