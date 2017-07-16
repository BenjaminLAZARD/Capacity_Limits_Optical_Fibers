function plot_proba_matrix( proba_matrix, Mrings, Nphases, Mrings_Rx, Nphases_Rx )

Colors = {'r','g','b','m','c','k','r','g','b','m','c','k','r','g','b','m','c','k'};

B=reshape(proba_matrix(1,:), [Nphases_Rx Mrings_Rx])';
C=reshape(proba_matrix(Nphases+1,:), [Nphases_Rx Mrings_Rx])';
D=reshape(proba_matrix(2*Nphases+1,:), [Nphases_Rx Mrings_Rx])';

B = sum(B');
C=sum(C');
D=sum(D');

%gaus = normpdf(1:Mrings_Rx, mean(B), sqrt(var(B)));

figure
plot(B, 'Color', Colors{1});
hold on
plot(C, 'Color', Colors{3});
hold on
plot(D, 'Color', Colors{4});
xlabel ('radii label');
ylabel ('probabilities');
title ('radii distribution of the first 3 input symbols');
legend( '1st symbol', '2nd symbol', '3rd symbol');

end

