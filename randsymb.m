function Y = randsymb( r, phi, N, M )
%% generates a matrix of random symbols from a constellation whose rings are at
%% vector r and phse at phi
		
     nr=numel(r);
     nphi=numel(phi);
     
     Ir=floor(rand(N,M)*nr)+1;
     Iphi=floor(rand(N,M)*nphi)+1;
     
     
     R=r(Ir);
     Phi=phi(Iphi);
      
     Y=R.*exp(j*Phi);
     
     if M==1
         Y=transpose(Y);
     end
     
      
      

end

