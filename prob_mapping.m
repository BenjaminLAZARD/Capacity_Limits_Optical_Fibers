function out = prob_mapping(points, dR, dphi);
	[I J] = size(points);
	out = zeros(I, J);
	for i = 1:I
		for j = 1:J
			a = points(i,j);
			r = round(abs(a)/dR); %no need to take precautions for overflow, since Max_ampl_Tx does that
			phi = round(angle(a)/dphi);
			P = phi*dphi;% true phase
			if P > 2*pi - dphi
				phi = 0;
            end
			out(i, j) = r*dR*exp(1i*P);
        end
    end
end