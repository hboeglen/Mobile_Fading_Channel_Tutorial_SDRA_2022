x=ones(25000,1);
rayleighchan = comm.RayleighChannel('SampleRate',20e3, 'MaximumDopplerShift',200, 'PathDelays', [0 200 600 1600 2400 5000]*1e-9, 'AveragePathGains', [-3 0 -2 -6 -8 -10],'Visualization','Impulse and frequency responses');
y = rayleighchan(x);
pause(5);
rayleighchan = comm.RayleighChannel('SampleRate',2e6, 'MaximumDopplerShift',200, 'PathDelays', [0 200 600 1600 2400 5000]*1e-9, 'AveragePathGains', [-3 0 -2 -6 -8 -10],'Visualization','Impulse and frequency responses', 'PathGainsOutputPort',true);
[y, gain] = rayleighchan(x);
pause(5);

t=0:1/2e6:(24999*(1/2e6));
f1=figure('position',[100 300 600 500]);
figure(f1),plot(t,20*log10(abs(gain(:,1))),'linewidth',3)
grid;
hold on;
for i=2:6
    figure(f1),plot(t,20*log10(abs(gain(:,i))), 'linewidth',3)
end
hold off
xlabel('Time (s)');
ylabel('Amplitude (dBV)');
f2=figure('position',[750 300 600 500]);
figure(f2),plot(real(gain(:,1)), imag(gain(:,1)), 'linewidth',3),hold on;
for i=2:6
    figure(f2),plot(real(gain(:,i)), imag(gain(:,i)), 'linewidth',3)
end
hold off;
grid;
xlabel('real (V)');
ylabel('imag (V)');
