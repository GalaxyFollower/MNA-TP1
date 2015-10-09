
function E = eigtest(bool,m,tol);
	if (bool == 1)
		E = eigAnalytic(m,0.51302,0.008,0.004,2,5.45);
	else
		E = eigNumeric(m,0.51302,0.008,0.004,2,5.45,tol);
	end		
end



