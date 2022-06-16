% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for simulating binary phase shift keyed transmission and
% reception and compare the simulated and theoretical bit error
% probability
% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 5 August 2007
% Adapted by Herve Boeglen for the Rayleigh channel 06/22
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
N = 2*10^6; % number of bits or symbols
rand('state',100); % initializing the rand() function
randn('state',200); % initializing the randn() function

% Transmitter
ip = rand(1,N)>0.5; % generating 0,1 with equal probability
s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 1 
n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white gaussian noise, 0dB variance 
Eb_N0_dB = [-2:2:10]; % multiple Eb/N0 values

for ii = 1:length(Eb_N0_dB)
   % Noise addition
   y = s + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise
   % receiver - hard decision decoding
   ipHat = real(y)>0;

   % counting the errors
   nErr(ii) = size(find([ip- ipHat]),2);

end

AWGNsimBer = nErr/N; % simulated ber

for ii = 1:length(Eb_N0_dB)
   % White Rayleigh process
   nr=(randn(1,N)+j*randn(1,N))*(1/sqrt(2));
   % Add noise
   y = s.*nr + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise
   y=y.*conj(nr);
   % receiver - hard decision decoding
   ipHat = real(y)>0;

   % counting the errors
   nErr(ii) = size(find([ip- ipHat]),2);

end

RayleighsimBer = nErr/N; % simulated ber
%theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
close all
f1=figure('position',[400 100 900 700]);
figure(f1)
semilogy(Eb_N0_dB,AWGNsimBer,'bx-','linewidth',3);
hold on
semilogy(Eb_N0_dB,RayleighsimBer,'mx-','linewidth',3);
axis([-2 10 10^-6 0.5])
grid on
legend('AWGN', 'Rayleigh');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK modulation');


