%**************************************************************************
%                Filtros FLG - M�todo de ventanas
%**************************************************************************
% Ejemplo: pasabajos con ventana rectangular 

clear all
close all

%% Par�metros de simulaci�n

nfft = 1024;  % tama�o de la FFT 
w = (0:nfft/2-1)/nfft*pi*2;

%% Especificaciones

ws = 0.2*pi;     % frecuencia de paso
wp = 0.32*pi;    % frecuencia de supresi�n
delta = 0.09;    % ripple 
DW = abs(ws-wp); % ancho de la banda de transici�n
wc = (ws+wp)/2;  % frecuencia de corte del filtro
M = ceil( 4*pi/DW );  % orden del filtro para ventana rectangular
n = 0:M;      % tiempo discreto

%% Filtro pasabajos

% Respuesta temporal del pasabajos con retardo de M/2.
hd = sinc(wc/pi*(n-M/2))*wc/pi;

% Ventaneo
win = window(@rectwin, M+1)';      % ventana rectangular
h = hd.*win;   % respuesta impulsiva del filtro ventaneado

% Respuesta en frecuencia del filtro
H = fft(h, nfft); H = H(1:nfft/2); 


%% Gr�ficos

% Respuesta impulsiva
stem(n, h, 'markerfacecolor','b')
xlabel('\Tiempo discreto')
ylabel('h[n]')
title('Respuesta impulsiva del filtro')
grid on

figure
% Respuesta en frecuencia
plot(w/pi,abs(H),'linewidth',2); 
title(['Filtro pasabajos, ventana rectangular. Orden M = ',num2str(M),', \omega_p = ',num2str(wp/pi),...
    '\pi, \omega_s = ',num2str(ws/pi),'\pi, \delta = ',num2str(delta),'.'])
ylabel('|H(e^{j\omega}|')
xlabel('\omega/\pi');
grid on

figure
% Diagrama de polos y ceros de la transferencia H(z)
b = h;  % coeficientes del numerador de H(z)
a = 1;  % coeficientes del denominador de H(z)
zplane(b,a);
title('Polos y ceros de H(z)')

figure
% retardos de fase
plot(w/pi, grpdelay(b,a,nfft/2) ,'LineWidth',2)
xlabel('\omega/\pi');
ylabel('Retardo [muestras]')
grid on

