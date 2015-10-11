function = eigPlot (en,ea) 
	plot(real(en),imag(en),'r.',real(ea),imag(ea),'b.')
	xlabel('Eje Real'); ylabel('Eje Imaginario'); legend('Eigenvalores numéricos','Eigenvalores Analíticos')
end
