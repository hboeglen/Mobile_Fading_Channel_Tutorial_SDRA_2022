pkg load signal;
clear all;
%==========================================================================
% PARAMETERS DEFINITION
%==========================================================================
fd=100; %Doppler frequency
fdmax=fd;
Ts = 10e-5; %Symbol time
sp=1000; %number of samples to be generated
      
%==========================================================================
% UNIT VARIANCE WHITE GAUSSIAN COMPLEX PROCESS GENERATION
%==========================================================================        
n=(randn(1,sp)+j*randn(1,sp))*(1/sqrt(2));

%==========================================================================
% JAKES FILTER DEFINITION IN THE FREQUENCY DOMAIN
%==========================================================================
fftlength = 1024; %FFT size
f = [0+eps:(fftlength/8)*(fd/fdmax)]/((fftlength/8)*(fd/fdmax));
for i=1:length(f)
    jpsd(i)=1/((1-f(i)^2)^.5);
    if (jpsd(i) > 100)
        jpsd(i)=100;
    end
end
psd=[jpsd(1:end) zeros(1,fftlength-2*(length(jpsd))) jpsd(end:-1:1)];
psdsqrt=sqrt(psd); %Square-root of the PSD ==> H(z)

%We get the IR by using an IFFT
filt=ifftshift(ifft(psdsqrt));
%Real coefficients
filt=real(filt).*hanning(1024)';%Hanning window
%Normalisation of the IR
filt=filt/sqrt(sum(filt.^2));

%==========================================================================
% FILTERING OF THE SAMPLES
%==========================================================================
path=fftfilt(filt,[n zeros(1,fftlength)]);
n = path(1+fftlength/2:end-fftlength/2);

%==========================================================================
% RESAMPLING
%==========================================================================
D_res = fd/2/0.0625; %Doppler filter resolution
%Symbol time = sampling time of the channel samples
O_res = 1/Ts
%Resampling factor calculation
m=round(O_res/D_res);
%Resampling
n_r=resample(n,m,1);
%==========================================================================
% RESULTS PRINTING
%==========================================================================
t=0:Ts:length(n_r)*Ts-Ts;
f1=figure('position',[100 300 600 500])
figure(f1),subplot(211),plot(t,20*log10(abs(n_r)))
axis([0 0.5 -40 10])
grid
xlabel('Time (s)')
ylabel('Amplitude (dBV)')
figure(f1),subplot(212),plot(t,angle(n_r))
axis([0 0.5 -4 4])
grid
xlabel('Time (s)')
ylabel('Phase (Rad)')
yy=abs(fft(n_r(1:end-10),65536*16));
yy=yy/max(yy);
f=(1/Ts)*[-65536*8:65536*8-1]/(65538*16);
f2=figure('position',[750 300 600 500])
figure(f2),plot(f,10*log10(fftshift(yy)))
grid
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
axis([-300 300 -30 0])



