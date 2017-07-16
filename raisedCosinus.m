function result = raisedCosinus(t, rho, T)
	result = ones(size(t));
	
	i = (t != 0);

	if (any (i))
		ratio = t(i)/T;
		a = sinc(ratio);
		b = cos(pi*rho*ratio);
		c = 1 - (2*rho*ratio).^2;
		result(i) = a .* b ./ c;
	endif
	
end