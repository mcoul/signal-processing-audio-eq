
config_m;


%% Parámetros de simulación
nfft = 1024;  % tamaño de la FFT 
w = (0:nfft/2-1)/nfft*pi*2; 
w = w/pi; % Normalizo
hEA = load('SEA.mat');
fs = 44100;
M = 80;

% Saco la data del hEA
h = hEA.h;


% Selector de ejs
bool_ej = [0 0 1];






if(bool_ej(1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ANÁLISIS DEL SEA %%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grafico la transformada

H=fft(h,nfft);
H=H(1:end/2);


plot(w,abs(H),'LineWidth',4)
xlabel('w')
ylabel('|H|')
grid minor
hold on

% Límites de interés

f_bottom = 20;
f_top = 16e3;

w_bottom = f_bottom/fs*2;
w_top = f_top/fs*2;

plot([w_top, w_top],[0,4],'color',[0 .5 0],'linewidth',3)
plot([w_bottom, w_bottom],[0,4],'color',[0 .5 0],'linewidth',4)


% Busco los extremos locales
[PKS, LOC, EXTRA]=findpeaks(abs(H'));
plot(LOC./nfft*2,PKS,'.r','MarkerSize',10)


[PKS, LOC2, EXTRA]=findpeaks(abs(1./H(1:f_top/fs*nfft-20)'));
plot(LOC2./nfft*2,abs(H(LOC2)),'r.','MarkerSize',10)

% Bandas de interés
LOC=sort(LOC);
LOC2=sort(LOC2);

w1 = (LOC(1)+LOC2(1))/2/(nfft/2);
w2 = (LOC2(2)+LOC(1))/2/(nfft/2);
w3 = (LOC(3)+LOC2(2))/2/(nfft/2);
w4 = (LOC(end)+LOC2(4))/2/(nfft/2);

wc = [w1 w2 w3 w4];

plot([w1 w1],[0 4],'y-.','Linewidth',3);
plot([w2 w2],[0 4],'y-.','Linewidth',3);
plot([w3 w3],[0 4],'y-.','Linewidth',3);
plot([w4 w4],[0 4],'y-.','Linewidth',3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DISEÑO DEL FILTRO MULTIBANDA %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aux=[H(LOC2(1)) H(LOC(1)) H(LOC2(2)) H(LOC(2)) H(LOC2(end))];
a = abs(1./aux);

% Teórico
%	wc = [0.05 0.12 0.203 0.456]
%	a = [2.33 0.284 3.93 0.569 1.2666];

% Copia teórica de HSEA
%	wc = [0.050 0.120 0.193 0.456]
%	a = aux;

%	a = [aux(1) aux(2)-0.6 aux(3) aux(4) aux(5)];
%	a = abs(1./a);

% Casi funciona con @hanning M=80 en vez de 40.
%	M = 80;
%	wc = [0.028 0.140 0.178 0.456]
%	a = [3.10 0.286 05.03 0.569 1.266];
%	tipo_ventana = @hanning;


% Funciona con @hamming M=90 en vez de 40.
	M = 90;
	wc = [0.028 0.140 0.178 0.456];
	a = [3.00 0.260 04.83 0.569 1.266];
	tipo_ventana = @hamming;


[heq, M] = multibanda(a,wc,M,tipo_ventana);
HEQ = fft(heq,nfft);
HEQ = HEQ(1:end/2);

figure
plot(w,abs(H),'LineWidth',4)
xlabel('w')
ylabel("Amplitud [\\si{\\dB}]")
grid minor
hold on
plot(w,abs(HEQ),'r','LineWidth',4)

plot([wc(1) wc(1)],[-00 05],'LineWidth',4,'y')
plot([wc(2) wc(2)],[-00 05],'LineWidth',4,'y')
plot([wc(3) wc(3)],[-00 05],'LineWidth',4,'y')
plot([wc(4) wc(4)],[-00 05],'LineWidth',4,'y')

%aux1=abs(fft(multibanda([a(1) 0],wc(1),M,tipo_ventana),nfft));
%plot(w,aux1(1:end/2),'color',[0 .5 0],'LineWidth',4)
%
%aux1=abs(fft(multibanda([0 a(2) 0],[wc(2) wc(1)],M,tipo_ventana),nfft));
%plot(w,aux1(1:end/2),'color',[0 .9 0],'LineWidth',4)


%% Gráfico de la rta del sistema ecualizado en dB

figure(3)
plot(w,mag2db(abs(HEQ.*H)),'LineWidth',4)
hold on
plot([w_top w_top],[-20 20],'LineWidth',4,'g')
plot(w,(w.*0)+2,'LineWidth',4,'r')
plot(w,(w.*0)-2,'LineWidth',4,'r')

lab_x = ["Frecuencia $\\frac{w}{\\pi}$"];
lab_y = ["Amplitud [\\si{\\dB}]"];
leyenda = ["Respuesta en frecuencia del sistema ecualizado por FIR FLG"];
leyenda = [leyenda; "Frecuencia máxima perceptible por el oido humano"];
loc = 'SouthWest';
AXIS = [0 1 -10 10];

set_graph('plot',[lab_x; lab_y],leyenda, loc, AXIS, 1);

plot([wc(1) wc(1)],[-20 20],'LineWidth',4,'y')
plot([wc(2) wc(2)],[-20 20],'LineWidth',4,'y')
plot([wc(3) wc(3)],[-20 20],'LineWidth',4,'y')
plot([wc(4) wc(4)],[-20 20],'LineWidth',4,'y')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(bool_ej(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Ejercicio 1.2 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 060;
dw = 0.025;
wc = [0.045 w2 w3 w4];
wc = [wc(1)-3e-3 w2+13e-3 w3-03e-3 w4];
w = [0 wc(1)-dw wc(1)+dw wc(2)-dw wc(2)+dw wc(3)-dw wc(3)+dw wc(4)-dw wc(4)+dw w_top-dw w_top+dw 1];
a = [2.33 0.284 3.93 0.569 1.2666 1];
a = [a(1) a(1) a(2) a(2) a(3) a(3) a(4) a(4) a(5) a(5) a(6) a(6)];
W = [35. 300 051 05 001.5 10];


heq2 = firls(M,w,a,W);
HEQ2 = fft(heq2,nfft)';
HEQ2 = HEQ2(1:end/2);

% Grafico la transformada
figure
w = (0:nfft/2-1)/nfft*pi*2; 
w = w/pi;

plot(w,abs(H),'LineWidth',4)
hold on
plot(w,abs(HEQ2),'r','LineWidth',4)
grid minor
plot([wc(1) wc(1)],[-00 05],'LineWidth',4,'y')
plot([wc(2) wc(2)],[-00 05],'LineWidth',4,'y')
plot([wc(3) wc(3)],[-00 05],'LineWidth',4,'y')
plot([wc(4) wc(4)],[-00 05],'LineWidth',4,'y')

%% Gráfico de la rta del sistema ecualizado en dB
figure
plot(w,mag2db(abs(HEQ2.*H)),'LineWidth',4)
hold on
plot([w_top w_top],[-20 20],'LineWidth',4,'g')
plot(w,(w.*0)+2,'LineWidth',4,'r')
plot(w,(w.*0)-2,'LineWidth',4,'r')

lab_x = ["Frecuencia $\\frac{w}{\\pi}$"];
lab_y = ["Amplitud [\\si{\\dB}]"];
leyenda = ["Respuesta en frecuencia del sistema ecualizado por cuadrados mínimos"];
leyenda = [leyenda; "Frecuencia máxima perceptible por el oido humano"];
loc = 'SouthWest';
AXIS = [0 1 -10 10];

set_graph('plot',[lab_x; lab_y],leyenda, loc, AXIS, 1);

plot([wc(1) wc(1)],[-20 20],'LineWidth',4,'y')
plot([wc(2) wc(2)],[-20 20],'LineWidth',4,'y')
plot([wc(3) wc(3)],[-20 20],'LineWidth',4,'y')
plot([wc(4) wc(4)],[-20 20],'LineWidth',4,'y')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(bool_ej(3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Ejercicio 2 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diagrama de polos y ceros
figure
b = h;	% Coeficientes del numerador de H(z)
a = 1;	% Coeficientes del denominador de H(z)
%zplane(b,a)


% Hallo los ceros de H.
ceros = roots(h);

% Destaco los que se encuentran afuera
aux = ceros.*(ceros>1);
aux2 = ceros-aux;

% Elimino los 0's redundantes
ceros_no_min =0;
ceros_min = 0;
for i=1:length(aux)
	if(aux(i))
		ceros_no_min = [ceros_no_min ; aux(i)];
	else
		ceros_min = [ceros_min ; aux2(i)];
	end
end
ceros_no_min = ceros_no_min(2:end);
ceros_min = ceros_min(2:end);

% Pongo los ceros dentro del círculo unitario
ceros_inv_conj = 1./conj(ceros_no_min);
fase_min = [ceros_min ; ceros_inv_conj];

h_fm = (poly(fase_min));

% Genero Hiir (compensador)
k = abs(h(1)*prod(ceros_no_min));
Hiir = 1./freqz(h_fm,1,nfft)/k;
Hiir = Hiir(1:end/2);

H = fft(h,nfft)';
H = H(1:end/2);

figure
%plot(mag2db(abs(H)),'LineWidth',4)
%hold on
%plot(mag2db(abs(Hiir)),'LineWidth',4,'g')
plot(w,mag2db(abs((Hiir).*H)),'LineWidth',4,'r')
grid minor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comentarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ej 20 de la guia 2 es parecido a lo que hay que haacer en 2.2.1
% 
