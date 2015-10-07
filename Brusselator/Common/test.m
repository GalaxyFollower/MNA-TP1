
function [E1,E2] = test(m,L,delta1,delta2,alpha,beta)
	E1 = eigBrusselatorJ(m,L,delta1,delta2,alpha,beta);
	E2 = eigBrusselatorJ4(m,L,delta1,delta2,alpha,beta);
endfunction

function [E1,E2]  =  test1 ()
	[E1,E2] = test(100,0.51302,0.008,0.004,2,5.45);
endfunction


function [E1,E2] = test2 ()
	[E1,E2] = test(1000,0.51302,0.008,0.004,2,5.45);
endfunction 


