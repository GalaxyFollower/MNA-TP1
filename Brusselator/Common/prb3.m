%% Tridiagonal QR algorithm
%% http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf

function E = eigBrusselatorJ3 (m,L,delta1,delta2,alpha,beta,EPSILON)
	A = brusselatorJMatrix(m,L,delta1,delta2,alpha,beta);
	i = 1; j = 2*m;
	ti = time();
	tw = 0;
	w = 0;
	qr_times = [];
	while (j>1)
		ti_while = time();
		A=A(1:j,1:j);
		ti_qr = time();
		[Q,R] = calculateQR (A, j);
		tf_qr = time() - ti_qr;
		qr_times(end + 1) = tf_qr;
		A = R * Q;
		if ( abs(A(j,j-1)) < EPSILON )
			E(i)=A(j,j);
			j=j-1;i=i+1;
		else
			Eaux = eig2p2(A(j-1:j,j-1:j));
			E(i) = Eaux(1); E(i+1) = Eaux(2); 
			i=i+2;j = j - 2;
		end
		tf_while = time() - ti_while;
		tw += tf_while;
		w++;
	end
	if(j==1)
		E(i)=A(1,1);
	end
	tf = time() - ti
	avg_while_time = tw / w
	qr_times
endfunction

function E = eigBrusselatorJ2 (m,L,delta1,delta2,alpha,beta)
	E = eigBrusselatorJ3 (m,L,delta1,delta2,alpha,beta,1e-10)
endfunction

function [Q,R] = calculateQR (A, m)
% each interaction of the k cicle lasts more than the previous one 
	for k=1:m
		Q(:,k) = A(:,k);
		for i=1:k-1
			Q(:,k) = Q(:,k) - A(:,k)*(transpose(Q(:,i)))*Q(:,i);
		end
		Q(:,k) = Q(:,k)/norm(Q(:,k));
	end
	R = transpose(Q)*A;
endfunction

function E = eig2p2 (A)
	p = [ 1 , A(1,1)+A(2,2) , A(1,1)*A(2,2) - A(1,2)+A(2,1) ];
	E = roots(p);
endfunction
