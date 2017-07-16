function I = mutual_inf( proba_matrix, Mrings, Nphases, Mrings_Rx, Nphases_Rx )

I=0;
for i=1:Mrings*Nphases
    for j=1:Mrings_Rx*Nphases_Rx
        P_Y=0;
        for k=1:Mrings*Nphases
            P_Y = P_Y + (1/(Mrings*Nphases))*proba_matrix(k,j);
        end
        if proba_matrix(i,j) ~= 0 &&  P_Y ~= 0
            I = I + (1/(Mrings*Nphases))*proba_matrix(i,j)*log2(proba_matrix(i,j)/P_Y);
        end
    end
end

