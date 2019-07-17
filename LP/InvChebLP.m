% Composition of Active and Passive Circuits 2019
% Inverse Chebyshev Low Pass Filter
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

%% Specifications
m=2;

f_p=1.1*(3+m)*1000
f_s=1.9*f_p
a_min=25+(max(1,a3)-5)*3/4
a_max=0.55+(max(1,a4)-5)/16
C_avail=0.1e-6
LowFreqGain=5

% f_p=4000
% f_s=8800
% a_min=22.75
% a_max=1
% C_avail=0.1e-6
% LowFreqGain=0

% f_p=1000/(2*pi)
% f_s=1400/(2*pi)
% a_min=18
% a_max=0.25
% C_avail=0.1e-6
% LowFreqGain=0

%% Frequency Regularization
w_p=2*pi*f_p
w_s=2*pi*f_s
W_p=w_p/w_s

%% Calculate System Rank
n=acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(1/W_p))
n=ceil(n)


epsilon=(10^(a_min/10)-1)^(-1/2)
a=1/n*asinh(1/epsilon)

%% Calculate Half Power Frequency
w_hp=1/cosh(1/n*acosh(1/epsilon))

%% Butterworth Angles
psi_k=zeros(1,n);

for k=1:n
    psi_k(k)=(pi*(2*k-1)/(2*n)-pi/2);
end

rad2deg(psi_k)

%% Calculate Poles
sigma_k=sinh(a)*cos(psi_k)
W_k=cosh(a)*sin(psi_k)
Wo_k=sqrt(sigma_k.^2+W_k.^2)
Q_k= 1./(2*cos(atan(W_k./sigma_k)))

%% Calculate Inverse Poles
wo_k=1./Wo_k;
p_k=1./(-sigma_k+i*W_k)'

%% Calculate Zeros
W_z=zeros(1,ceil(n/2));
for k=1:ceil(n/2)
    W_z(k)=sec((2*k-1)*pi/(2*n));
end
if(ceil(n/2)~=n/2)
    W_z(end+1)=Inf;
end
wo_k'
W_z

%% Low Pass Notch

N=cell(ceil(n/2),1);
D=cell(ceil(n/2),1);
T=cell(ceil(n/2),1);

for i=1:ceil(n/2)
    if imag(p_k(i))~=0
        omega_z(i)=W_z(i)/wo_k(i)
        R_1old(i)=1;
        R_4old(i)=1;
        C_old(i)=1/(2*Q_k(i))
        R_2old(i)=4*Q_k(i)^2
        R_5old(i)=4*Q_k(i)^2/(omega_z(i)^2-1)
        R_3old(i)=omega_z(i)^2/(2*Q_k(i)^2)
        kH(i)=1/(1+omega_z(i)^2/(2*Q_k(i)^2))
        kL(i)=kH(i)*omega_z(i)^2
        k_f(i)=w_s*wo_k(i)
        k_m(i)=C_old(i)/(k_f(i)*C_avail)
        
        C(i)=C_avail
        R_1(i)=R_1old(i)*k_m(i)
        R_2(i)=R_2old(i)*k_m(i)
        R_3(i)=R_3old(i)*k_m(i)
        R_4(i)=R_4old(i)*k_m(i)
        R_5(i)=R_5old(i)*k_m(i)
        
        N{i}=kH(i)*[1 , ((kH(i)-1)/kH(i)*1/(R_1(i)*C(i)) +2/(R_2(i)*C(i))+2/(R_5(i)*C(i)) ) , (1/(R_1(i)*R_5(i)*C(i)^2 )+1/(R_1(i)*R_2(i)*C(i)^2 ))];
        D{i}=[1, 2/(R_2(i)*C(i)) ,1/(R_1(i)*R_2(i)*C(i)^2 )];
        T{i}=tf(N{i},D{i});
    elseif imag(p_k(i))==0
        C_old(i)=1
        R_old=1/wo_k(i)
        k_f(i)=w_s
        k_m(i)=C_old(i)/(k_f(i)*C_avail)
        C(i)=C_avail
        R=k_m(i)*R_old
        N{i}=1/(R*C(i));
        D{i}=[1, 1/(R*C(i))];
        T{i}=tf(N{i},D{i});
    end
end

for i=1:ceil(n/2)
    plot_transfer_function(T{i},[f_s,f_p])
    name = ['pics/T',num2str(i),'.png'];
    saveas(gcf,name);
end

k=-10^(LowFreqGain/20)/prod(kL)
r1=1e4
r2=-k*r1
T_LP=tf(k,1);
for i=1:ceil(n/2)
    T_LP=T{i}*T_LP;
end

plot_transfer_function(T_LP,[f_s,f_p])
saveas(gcf,'pics/T_LP.png');

figure('Position', get(0, 'Screensize'));
for i=1:ceil(n/2)
    bodemag(T{i});
    hold on;
end
bodemag(T_LP);
grid on;
legend('Unit 1','Unit 2','T_{LP}')
saveas(gcf,'pics/bodeALL.png');

plot_transfer_function(inv(T_LP),[f_s,f_p])
saveas(gcf,'pics/invLP.png');

plot_transfer_function(1/10^(LowFreqGain/20)*T_LP,[f_s,f_p])
saveas(gcf,'pics/checkSpecs.png');

%% Spectrum Calculation of Input
input= @(t) 0.5*square(2*pi*2000*t,40)+0.5;
spectrum(T_LP,input)

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