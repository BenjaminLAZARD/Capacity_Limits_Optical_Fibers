function y = sinc(x)
	
	y=sin(pi*x)./(pi*x);
	
	I=find(x==0);
	
	y(I)=1;
	