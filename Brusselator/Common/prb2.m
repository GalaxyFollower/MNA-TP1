function E = eigBrusselatorJ (m,L,delta1,delta2,alpha,beta)
	h = 1/(m+1); 
	tau1  = delta1/(h*L)^2; 
	tau2  = delta2/(h*L)^2;
    	for j=1:m,
       		eigofT(j) = -2*(1- cos(pi*j*h) );  % eigenvalues of T
    	end;
    	for j=1:m,
      		coeff(1) = 1;
       		coeff(2) = alpha^2 - (beta - 1) - (tau1+tau2)*eigofT(j);
       		coeff(3) = beta*alpha^2 + tau1*tau2*eigofT(j)^2 + tau2*(beta-1)*eigofT(j) - alpha^2*tau1*eigofT(j) - alpha^2*(beta-1);
       		d = roots(coeff);
       		Aeig(j) = d(1);  Aeig(m+j) = d(2); % eigenvalues of A
    	end;
endfunction
