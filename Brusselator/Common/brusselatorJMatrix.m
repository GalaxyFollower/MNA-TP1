function A = brusselatorJMatrix (m,L,delta1,delta2,alpha,bet)
%%  PARAMETERS
%%  ----------
%%  m: 1/2 order of the matrix
%%  L: bifurcation parameter
%%  delta1: diffusion coefficient for x
%%  delta2: diffusion coefficient for y
%%  alpha: constant in reaction term for x
%%  bet: constant in reaction term for y
	h = 1/(m+1);
	aux = (h*L)^2;
	tau1 = delta1 / aux;
	tau2 = delta2 / aux;
	I = eye(m);
	T = tridiag(m,1,-2,1);
	M11 = tau1 * T + (bet-1)*I;
	M12 = (alpha^2) * I;
	M21 = -bet * I;
	M22 = tau2 * T - (alpha^2) * I;
	A = [ M11 , M12 ; M21 , M22 ];
end


function T = tridiag (dim,a,b,c)
	T = zeros(dim,dim);
	for k=1:1:dim
		T(k,k)=b;
		if ( k < dim )
			T(k,k+1)=c;
			T(k+1,k)=a;
		end
	end
end

