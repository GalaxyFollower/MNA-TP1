function en1 = eigPlot (en1,en2,en3,en4,ea) 
	plot(real(en1),imag(en1),'r.',real(en2),imag(en2),'g.',real(en3),imag(en3),'c.',real(en4),imag(en4),'m.',real(ea),imag(ea),'b.')
	xlabel('Eje Real'); ylabel('Eje Imaginario'); legend('Numericos con tolerancia 1e-4','Numericos con tolerancia 1e-6','Numericos con tolerancia 1e-8','Numericos con tolerancia 1e-10','Analiticos')
end
