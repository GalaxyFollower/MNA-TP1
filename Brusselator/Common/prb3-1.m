                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        function E = eigBrusselatorJ4 (m,L,delta1,delta2,alpha,beta)
	tol = 1E-10;
	A = brusselatorJMatrix(m,L,delta1,delta2,alpha,beta);
	A = hessenberg(A); % A(hess) = Q'*A*Q
	n = (size(A))(1);
	i=1; tol = 1E-10;
        while (n > 2)
		A(1:n,1:n);
                [Q,R]=calculateQR1(A);
		A = R*Q;
		if ( abs(A(n,n-1)) > tol*(abs(A(n-1,n-1))+abs(A(n,n))) )  %single shift
			E(i) = A(n,n);
			i=i+1; n = n-1;		
		elseif ( abs(A(n,n-1)) > tol*(abs(A(n-1,n-1))+abs(A(n,n))) )  %double shift
			Eaux = eig2p2 (A(n-1:n,n-1:n));
                        E(i) = Eaux(1);
			 E(i+1) = Eaux(2);
                        i=i+2;n = n - 2;
		endif	
        endwhile
	if (n==2)
		Eaux = eig2p2 (A);
                E(i) = Eaux(1); E(i+1) = Eaux(2);
        elseif (n==1)
		E(i) = A(1,1);
	endif
endfunction



function H = hessenberg(A) %IMPLEMENTAR
	H = hess(A);
endfunction


function [Q,R] = calculateQR1 (A)
% NO FUNCIONA
	n = (size(A))(1);
        for k=1:n
                Q(:,k) = A(:,k);
                for i=1:k-1
			aux = transpose(A(:,k))*Q(:,i);
                        Q(:,k) = Q(:,k) - aux*Q(:,i)
                end
                Q(:,k) = Q(:,k)/norm(Q(:,k));
        endfor
        R = transpose(Q)*A;
endfunction




function E = eig2p2 (A)
        p = [ 1 , -A(1,1)-A(2,2) , A(1,1)*A(2,2) - A(1,2)*A(2,1) ];
        E = roots(p);
endfunction
