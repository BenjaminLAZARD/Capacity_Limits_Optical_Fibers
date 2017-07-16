function [ rings, phases, P, s , svec] = constel( Mrings, Nphases, Max_ampl )

dring = Max_ampl/Mrings; %  uniformly distributed rings
dphases = 2*pi/Nphases; % uniformlfy distributed phases
%%% Generate the Tx Ring constellation
rings = dring:dring:Max_ampl;
phases = 0:dphases:2*pi-dphases; % uniformly distributed phases


s = zeros(Mrings, Nphases);
for i = 1:Mrings
	s(i,:) = rings(i)*exp(1i*phases);
	%figure(1)
	%plot(s(i,:), 'o');
end

P=mean(rings.^2);
svec=reshape(s',[1,Mrings*Nphases]);
end



