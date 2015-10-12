function E = eigAnalytic (m,L,delta1,delta2,alpha,bet)
%%  PARAMETERS
%%  ----------
%%  m: 1/2 order of the matrix
%%  L: bifurcation parameter
%%  delta1: diffusion coefficient for x
%%  delta2: diffusion coefficient for y
%%  alpha: constant in reaction term for x
%%  bet: constant in reaction term for y
	h = 1/(m+1); 
	tau1  = delta1/(h*L)^2; 
	tau2  = delta2/(h*L)^2;
    	for k=1:m,
       		eigofT(k) = -2*(1- cos(pi*k*h) );  % eigenvalues of T
    	end
    	for k=1:m
      		coeff(1) = 1;
       		coeff(2) = alpha^2 - (bet - 1) - (tau1+tau2)*eigofT(k);
       		coeff(3) = bet*alpha^2 + tau1*tau2*eigofT(k)^2 + tau2*(bet-1)*eigofT(k) - alpha^2*tau1*eigofT(k) - alpha^2*(bet-1);
       		d = roots(coeff);
       		E(k) = d(1);  E(m+k) = d(2); % eigenvalues of A
    	end
end
