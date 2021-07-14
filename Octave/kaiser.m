

function [w, M] = kaiser (A, dw)

	% Calculo de beta
	if A > 50
		beta = 0.1102*(A-8.7);
        elseif A >= 21
		beta = 0.5842*(A-21)^(0.4) + 0.07886*(A-21);
        else
		beta = 0;
	end

	% Calculo de M
	M = ceil((A-8)/(2.285*dw*pi));
    
	% Calculo de w(n)
	n = 0:M;	% Genero un vector de 0 a M
	w = n.*0;
	alpha = (M)/2;
    
        w = besseli(0,beta*sqrt(1-((n-1-alpha)/alpha).^2))./besseli(0,beta);
	w = w';
 

endfunction
