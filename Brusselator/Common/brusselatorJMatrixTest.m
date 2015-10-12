function b = brusselatorJMatrixTest()

	A200 = mmread('bwm200.mtx');
	A200 = full(A200);
	B200 = brusselatorJMatrix(100,0.51302,0.008,0.004,2,5.45);
	aux = all(A200-B200);
	b='success';
	for k=1:length(aux)
		if (aux(k)!=0)
			b='failed bwm200.mtx';
		end
	end
	if strcmp(b,'success')
		A2000 = mmread('bwm2000.mtx');
        	A2000 = full(A2000);
        	B2000 = brusselatorJMatrix(1000,0.51302,0.008,0.004,2,5.45);
        	aux = all(A2000-B2000);
        	for k=1:length(aux)
                	if (aux(k)!=0)
                        	b='failed bwm2000.mtx';
                	end
        	end
	end
end
