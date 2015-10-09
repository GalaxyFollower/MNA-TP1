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

% ---------------- HOUSE HOLDER QR FUNCTIONS ---------------- %

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

function [Q,R] = calculateHQR2(A)
    n = length(A); 
    m= max(abs(A)); 
    A = A/m; 
    q = A(2:n)'* A(2:n); 
    Q(1)=1; 
    Q(2:n)=x(2:n); 
    if q == 0 
        R=0; 
    else 
        a = (A(1)^2+q)^(1/2); 
        if x(1) <= 0 
            Q(1) = x(1)-a; 
        else 
            Q(1)=-q/(x(1)+a); 
        end 
        R = 2*Q(1)^2 / (q+Q(1)^2); 
        Q = Q/Q(1); 
    end 
    Q = Q'; 
end

function [Q,R] = calculateHQR3(A)
    [m, n] = size(A);
    R = A;
    for k = 1:n,
      x = R(k:m,k);
      e = zeros(length(x),1); e(1) = 1;
      Q = sign(x(1))*norm(x)*e + x;
      Q = Q. /norm(Q);
      R(k:m, k:n) = R(k:m, k:n) -2*Q*Q’*R(k:m, k:n);
      Q(k:m,k) = Q;
end

% ---------------- HESSENBERG FUNCTIONS ---------------- %

function H = hessenberg(A) 
    H = hess(A);
end

function [L,H] = hessenberg2(A) 
% si la matriz no es cuadrada explota    
    [m,n] = size(A); 
    L = zeros(m,n);
    H = A;
    for j = 1:m-2
        x = H(j+1:m,j);
        x(1) = x(1) + sign(x(1)) * norm(x);
        n = norm(x);
        if n > 0
            u = x/norm(x);
            H(j+1:m,j:m) = H(j+1:m,j:m) - 2*u*( u'*H(j+1:m,j:m) );
            H(1:m,j+1:m) = H(1:m,j+1:m) - 2*( H(1:m,j+1:m)*u )*u';
        else
            u = x;
        end
    L(j+1:m,j) = u;
  end
end

function [ H, Q ] = hessenberg3(A)
    [m, n] = size(A);
    if issparse(A)
        I = speye(size(A));
        v = sparse(zeros(m, 1));
    else
        I = eye(size(A));
        v = zeros(m, 1);
    end

    PP = I;

    for k=1:m-2
        alpha = -sign(A(k+1, k)) * sqrt(sum(A(k+1:end, k).^2));
        r = sqrt((alpha^2-A(k+1, k)*alpha)/2);
        v(1:k) = zeros(k, 1);
        v(k+1) = (A(k+1, k)-alpha)/(2*r);
        v(k+2:m) = A(k+2:m, k)./(2*r);
        P = I - 2*v*v';
        PP = PP * P;
        A = P*A*P;
    end

    H = triu(A, -1);
    Q = PP;
end

% ---------------- END HESSENBERG FUNCTIONS ---------------- %

function E = eig2p2 (A)
    p = [ 1 , -A(1,1)-A(2,2) , A(1,1)*A(2,2) - A(1,2)*A(2,1) ];
    E = roots(p);
end
