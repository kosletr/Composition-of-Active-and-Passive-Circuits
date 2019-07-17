% Composition of Active and Passive Circuits 2019
% Chebyshev High Pass Filter
% Letros Konstantinos 8851

%% Clear
clear
clc

format long g

%% My Info

a1=8;
a2=8;
a3=5;
a4=1;

m=0;

%% Specifications
f_p=(3+m)*1000
f_s=f_p/1.8
a_min=25+a3*4/9
a_max=0.5+a4*0.25/9
C_avail=0.1e-6
LowFreqGain=0

% f_p=25120/(2*pi)
% f_s=13941.6/(2*pi)
% a_min=21.11
% a_max=0.63
% C_avail=1e-6
% LowFreqGain=0

% f_p=7500
% f_s=5000
% a_min=22
% a_max=0.5
% C_avail=0.1e-6
% LowFreqGain=0

%% Calculate System's Rank
w_p=2*pi*f_p
w_s=2*pi*f_s
W_s=w_p/w_s

n=acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(W_s))
n=ceil(n)

%% Calculate Half Power Frequency
epsilon=(10^(a_max/10)-1)^(1/2)
a=1/n*asinh(1/epsilon)
W_hp=cosh(1/n*acosh(1/epsilon))

%% Calculate Butterworth Angles
psi_k=zeros(1,n);
for k=1:n
    psi_k(k)=(pi*(2*k-1)/(2*n)-pi/2);
end

rad2deg(psi_k)'

%% Calculate Poles
Sigma_k=(sinh(a)*cos(psi_k))';
W_k=(cosh(a)*sin(psi_k))';
wo_k=sqrt(Sigma_k.^2+W_k.^2)
Q_k=wo_k./(2*Sigma_k)

p_k=(-Sigma_k+i*W_k)

%% Calculate Inverse Poles
w_hp=w_p/W_hp

for i=1:n
    if(imag(p_k(i))<=0)
        w_k(i)=w_p/abs(p_k(i));
        Q(i)=Q_k(i);
    end
end
w_k'
Q'

%% High Pass Sallen-Key (2nd Strategy)

N=cell(ceil(n/2),1);
D=cell(ceil(n/2),1);
T=cell(ceil(n/2),1);

for i=1:ceil(n/2)
    if imag(p_k(i))~=0
        C_old(i)=1
        k(i)=1
        R_1old(i)=1/(2*Q(i))
        R_2old(i)=2*Q(i)
        k_f(i)=w_k(i)
        k_m(i)=C_old(i)/(k_f(i)*C_avail)
        
        C(i)=C_avail
        R_1(i)=R_1old(i)*k_m(i)
        R_2(i)=R_2old(i)*k_m(i)
        
        N{i}=[1,0,0];
        D{i}=[1, w_k(i)/Q(i) ,w_k(i)^2];
        T{i}=tf(N{i},D{i});
    elseif imag(p_k(i))==0
        C_old(i)=1
        R_old=1
        
        k_f(i)=w_k(i)
        k_m(i)=C_old(i)/(k_f(i)*C_avail)
        C(i)=C_avail
        R=k_m(i)*R_old
        
        N{i}=[1,0];
        D{i}=[1, 1/(R*C(i))];
        T{i}=tf(N{i},D{i});
    end
end

for i=1:ceil(n/2)
    plot_transfer_function(T{i},[f_s,f_p])
    name = ['pics/T',num2str(i),'.png'];
    saveas(gcf,name);
end

K=-10^(LowFreqGain/20)
r1=1e4
r2=-K*r1
T_HP=tf(K,1);
for i=1:ceil(n/2)
    T_HP=T{i}*T_HP;
end

plot_transfer_function(T_HP,[f_s,f_p])
saveas(gcf,'pics/T_HP.png');

figure('Position', get(0, 'Screensize'));
for i=1:ceil(n/2)
    bodemag(T{i});
    hold on;
end
bodemag(T_HP);
grid on;
legend('Unit 1','Unit 2','Unit 3','T_{HP}')
saveas(gcf,'pics/bodeALL.png');

plot_transfer_function(inv(T_HP),[f_s,f_p])
saveas(gcf,'pics/invHP.png');

plot_transfer_function(1/10^(LowFreqGain/20)*T_HP,[f_s,f_p])
saveas(gcf,'pics/checkSpecs.png');

%% Spectrum Calculation of Input
f1=0.4*w_s
f2=0.9*w_s
f3=1.4*w_p
f4=2.4*w_p
f5=4.5*w_p

input= @(t) cos(f1*t)+0.5*cos(f2*t)+cos(f3*t)+0.7*cos(f4*t)+0.5*cos(f5*t);
spectrum(T_HP,input)

function spectrum(sys,input)
T = 1e-2;
F_s = 1e6;
dt = 1/F_s;
t = 0:dt:T-dt;

input=input(t);
figure('Position', get(0, 'Screensize'));
plot(t,input);
axis([0 inf -0.5 1.5]);
title('Input Signal');
xlabel('Time in sec');
saveas(gcf,'pics/input.png');

y=lsim(sys,input,t);
figure('Position', get(0, 'Screensize'));
plot(t,y);
title('Output Signal');
xlabel('Time in sec');
saveas(gcf,'pics/output.png');

figure('Position', get(0, 'Screensize'));
plot(t,input);
hold on;
plot(t,y);
hold off;
title('Input and Output sigmals');
xlabel('Time in sec');
saveas(gcf,'pics/inputOutput.png');

Xf=fft(input);
L=length(input);

P2 = abs(Xf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = F_s*(0:(L/2))/L;


figure('Position', get(0, 'Screensize'));
plot(f,P1);
axis([0.01 20000 0 inf]);
title('FFT of Input Signal');
xlabel('Time in sec');
Yf=fft(y);
L=length(y);
saveas(gcf,'pics/FFTinput.png');

P2 = abs(Yf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = F_s*(0:(L/2))/L;
figure('Position', get(0, 'Screensize'));
plot(f,P1);
axis([0.01 20000 0 inf]);
title('FFT of Output Signal');
saveas(gcf,'pics/FFToutput.png');

end