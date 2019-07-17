function spectrum(sys,input)
T = (1/100);
Fs = 1e6;
dt = 1/Fs;
t = 0:dt:T-dt;

input=input(t);
figure;
plot(t,input);
axis([0 inf -0.5 1.5]);
title('Input signal');

y=lsim(sys,input,t);
figure;
plot(t,y);
title('Output signal');

figure;
plot(t,input);
hold on;
plot(t,y);
hold off;
title('Input and Output sigmals');



Xf=fft(input);
L=length(input);

P2 = abs(Xf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('FFT of Input signal');
Yf=fft(y);
L=length(y);

P2 = abs(Yf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('FFT of Output signal');

end