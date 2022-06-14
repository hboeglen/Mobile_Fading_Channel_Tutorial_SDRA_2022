pkg load signal;
clear all;
%==========================================================================
% PARAMETERS DEFINITION
%==========================================================================
fd=100; %Frequence Doppler
fdmax=fd;
Ts = 10e-5; %Temps symbole
sp=1000; %Longueur du vecteur representant les echantillons canal
        %avant reechantillonnage

%==========================================================================
% UNIT VARIANCE WHITE GAUSSIAN COMPLEX PROCESS GENERATION
%==========================================================================        
n=(randn(1,sp)+j*randn(1,sp))*(1/sqrt(2));

%==========================================================================
% JAKES FILTER DEFINITION IN THE FREQUENCY DOMAIN
%==========================================================================
fftlength = 1024; %longueur de la FFT
f = [0+eps:(fftlength/8)*(fd/fdmax)]/((fftlength/8)*(fd/fdmax));
for i=1:length(f)
    jpsd(i)=1/((1-f(i)^2)^.5);
    if (jpsd(i) > 100)
        jpsd(i)=100;
    end
end
psd=[jpsd(1:end) zeros(1,fftlength-2*(length(jpsd))) jpsd(end:-1:1)];
psdsqrt=sqrt(psd); %On prend la racine de la DSP ==> H(z)

%On obtient la RI du filtre en temporel par une IFFT
filt=ifftshift(ifft(psdsqrt));
%On veut un filtre de coefficients reels
filt=real(filt).*hanning(1024)';%Fenetrage de Hanning
%Normalisation de la RI
filt=filt/sqrt(sum(filt.^2));

%==========================================================================
% FILTERING OF THE SAMPLES
%==========================================================================
path=fftfilt(filt,[n zeros(1,fftlength)]);
n = path(1+fftlength/2:end-fftlength/2);

%==========================================================================
% RESAMPLING
%==========================================================================
D_res = fd/2/0.0625; %Resolution du filtre Doppler
%Temps symbole = temps d'�chantillonnage des echantillons canal
O_res = 1/Ts
%Calcul du facteur de surechantillonnage
m=round(O_res/D_res);
%Reechantillonnage
n_r=resample(n,m,1);
%==========================================================================
% RESULTS PRINTING
%==========================================================================
t=0:Ts:length(n_r)*Ts-Ts;
figure(1),subplot(211),plot(t,20*log10(abs(n_r)))
axis([0 0.5 -40 10])
grid
xlabel('Time (s)')
ylabel('Amplitude (dBV)')
figure(1),subplot(212),plot(t,angle(n_r))
axis([0 0.5 -4 4])
grid
xlabel('Time (s)')
ylabel('Phase (Rad)')
yy=abs(fft(n_r(1:end-10),65536*16));
yy=yy/max(yy);
f=(1/Ts)*[-65536*8:65536*8-1]/(65538*16);
figure(2),plot(f,10*log10(fftshift(yy)))
grid
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
axis([-300 300 -30 0])



