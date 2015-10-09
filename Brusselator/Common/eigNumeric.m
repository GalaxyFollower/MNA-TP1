function E = eigNumeric (m,L,delta1,delta2,alpha,bet,tol)
	A = brusselatorJMatrix(m,L,delta1,delta2,alpha,bet);
	A = hessenberg(A); % A(hess) = Q'*A*Q
	n = length(A);
	i=1;
    while (n > 2)
		A(1:n,1:n);
        [Q,R]=calculateHQR(A);
		A = R*Q;
		if ( abs(A(n,n-1)) < tol*(abs(A(n-1,n-1))+abs(A(n,n))) )  %single shift
			E(i) = A(n,n);
			i=i+1; n = n-1;		
		elseif ( abs(A(n-1,n-2)) < tol*(abs(A(n-1,n-1))+abs(A(n-2,n-2))) )  %double shift
			Eaux = eig (A(n-1:n,n-1:n));
            E(i) = Eaux(1); E(i+1) = Eaux(2);
            i=i+2;n = n - 2;
		end
    end
	if (n==2)
		Eaux = eig2p2 (A);
        E(i) = Eaux(1); E(i+1) = Eaux(2);
    elseif (n==1)
		E(i) = A(1,1);
	end
end

function [Q,R] = calculateHQR (A)
	%Es de Lucas, implementar nosotros, pero sabemos que anda
	m = length(A);
    Q = eye(m);
	R=A;
    for k=1:m-1
        u = R(k:end,k);
        
        u(1) = u(1) + sign(u(1))*norm(u);
        u = u/norm(u);

		R(k:end,k:end)=R(k:end,k:end)-(2*u)*(u'*R(k:end,k:end));
        Q(:,k:end)=Q(:,k:end)-(Q(:,k:end)*(2*u))*u';
        
        %Para transformar los coef de R en positivos.
        if(R(k,k)<0)
            R(k,k:end)=-R(k,k:end);
            Q(:,k)=-Q(:,k);
        end

    end

    for i=1:m-1
        for j=1:i
            R(i+1,j)=0;
        end
    end
end

function H = hessenberg(A) %IMPLEMENTAR
	H = hess(A);
end


function [Q,R] = calculateQR (A)
	% Ya no lo usamos
	n = (size(A))(1);
    for k=1:n
        Q(:,k) = A(:,k);
        for i=1:k-1
			aux = transpose(A(:,k))*Q(:,i);
            Q(:,k) = Q(:,k) - aux*Q(:,i)
        end
        Q(:,k) = Q(:,k)/norm(Q(:,k));
    end
    R = transpose(Q)*A;
end




function E = eig2p2 (A)
    p = [ 1 , -A(1,1)-A(2,2) , A(1,1)*A(2,2) - A(1,2)*A(2,1) ];
    E = roots(p);
end
