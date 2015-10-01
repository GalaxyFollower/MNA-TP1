%% Tridiagonal QR algorithm
%% http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf

function E = eigBrusselatorJ2 (m,L,delta1,delta2,alpha,beta,EPSILON)
	A = brusselatorJMatrix(m,L,delta1,delta2,alpha,beta);
	i = 1; j = m;
	while (j>1)
		A=A(1:j,1:j);
		[Q,R] = calculateQR (A);
		A = R * Q;
		if ( abs(A(j,j-1)) < EPSILON )
			E(i)=A(m-i,m-i);
			j=j-1;i=i+1;
		else
			Eaux = eig2x2(A(j-2:j,j-2:j);
			E(i) = Eaux(1); E(i+1) = Eaux(2); 
			i=i+2;j = j - 2;
		end
	end
	if(j=1)
		E(i)=A(1,1);
	end
endfunction



function E = eigBrusselatorJ2 (m,L,delta1,delta2,alpha,beta)
	E = eigBrusselatorJ2 (m,L,delta1,delta2,alpha,beta,10^-10)
endfunction



function [Q,R] = calculateQR (A) 
	for k=1:(2*m)
		Q(:,k) = A(:,k);
		for i=1:k-1
			Q(:,k) = Q(:,k) - A(:,k)*((Q(:,i))')*Q(:,i);
		end
		Q(:,k) = Q(:,k)/norm(Q(:,k));
	end
	R = Q'*A;
endfunction




function E = eig2x2 (A)
	p = [ 1 , A(1,1)+A(2,2) , A(1,1)*A(2,2) - A(1,2)+A(2,1) ];
	E = roots(p);
endfunction
