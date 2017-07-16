function matrix = prob_mat( prob, Mrings, Nphases, Mrings_Rx, Nphases_Rx, dring_Rx, dphases_Rx )

[M, N]= size(prob);
matrix = zeros(Mrings*Nphases, Mrings_Rx*Nphases_Rx);

for i=1:Mrings
    for k=1:N
        point = prob(i, k);
        ang = angle(point);
        if ang < 0
            ang = -ang;
        end
        if round(abs(point)/dring_Rx) == 0
            distance = 1;
        else
            distance = round(abs(point)/dring_Rx);
        end
        
        number = (distance-1)*Nphases_Rx + 1 + round((ang/2*pi)*Nphases_Rx);
        if number >=Nphases_Rx*Mrings_Rx
            number = 65536;
        end    
        
        matrix(Nphases*(i-1)+1, number) = matrix(Nphases*(i-1)+1, number) + 1/N;
    end
    for j=2:Nphases
        matrix((i-1)*Mrings + j, :) = circshift(matrix(Nphases*(i-1)+1,:), [1 Nphases_Rx/Nphases*(j-1)]);
    end
end