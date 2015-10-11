function E = eigNumeric (m,L,delta1,delta2,alpha,bet,tol)
%%  PARAMETERS
%%  ----------
%%  m: 1/2 order of the matrix
%%  L: bifurcation parameter
%%  delta1: diffusion coefficient for x
%%  delta2: diffusion coefficient for y
%%  alpha: constant in reaction term for x
%%  bet: constant in reaction term for y
	A = brusselatorJMatrix(m,L,delta1,delta2,alpha,bet);
	A = hessenberg(A);
	n = length(A);
	i=1;
	while (n > 2)
		[Q,R]=calculateHQR(A);
		A = R*Q;
		if ( abs(A(n,n-1)) < tol*(abs(A(n-1,n-1))+abs(A(n,n))) )  %single shift
			E(i) = A(n,n);
			i=i+1; n = n-1;	
			A = A(1:n,1:n);	
		elseif ( abs(A(n-1,n-2)) < tol*(abs(A(n-1,n-1))+abs(A(n-2,n-2))) )  %double shift
			Eaux = eig (A(n-1:n,n-1:n));
			E(i) = Eaux(1); E(i+1) = Eaux(2);
			i=i+2;n = n - 2;
			A = A(1:n,1:n);
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
	[m,n] = size(A);
	Q = eye(m);
	R = A;
	for j = 1:n
		normx = norm(R(j:end,j));
		s = -sign(R(j,j));
		u1 = R(j,j) - s*normx;
		w = R(j:end,j)/u1;
		w(1) = 1;
		tau = -s*u1/normx;
		R(j:end,:) = R(j:end,:)-(tau*w)*(w'*R(j:end,:));
		Q(:,j:end) = Q(:,j:end)-(Q(:,j:end)*w)*(tau*w)';
	end
end

function H = hessenberg(A)
	[m,n] = size(A); 
	L = zeros(m,n);
	H = A;
	for j = 1:m-2
		x = H(j+1:m,j);
		x(1) = x(1) + sign(x(1)) * norm(x);
		n = norm(x);
		if n > 0
			u = x/norm(x);
			H(j+1:m,j:m) = H(j+1:m,j:m) - 2*u*(u'*H(j+1:m,j:m));
			H(1:m,j+1:m) = H(1:m,j+1:m) - 2*( H(1:m,j+1:m)*u )*u';
		else
			u = x;
		end
        %L(j+1:m,j) = u;
        end
end


function [Q,R] = calculateQR (A)
	% Ya no lo usamos
	n = length(A);
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
