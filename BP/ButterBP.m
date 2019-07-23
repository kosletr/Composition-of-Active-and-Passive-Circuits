% Composition of Active and Passive Circuits 2019
% Butterworth Band Pass Filter
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

f_0=1200
f_1=650+25*a4
f_2=f_0^2/f_1
D=2.8*(f_0^2-f_1^2)/f_1;
f_3=(-D+sqrt(D^2+4*f_0^2))/2
f_4=f_0^2/f_3
a_min=24+a3*5/9
a_max=0.5+a4/36
C_avail=0.1e-6
LowFreqGain=10;

%  f_1=700
%  f_2=1428.571
%  f_3=408.426
%  f_4=2448.424
% a_min=25
% a_max=0.722
% C_avail=0.1e-6
% LowFreqGain=0

%  f_1=750/(2*pi)
%  f_2=1500/(2*pi)
%  f_3=300/(2*pi)
%  f_4=3750/(2*pi)
% a_min=22
% a_max=0.5
% C_avail=0.1e-6
% LowFreqGain=0

%% Calculate bandwidth
w_1=2*pi*f_1
w_2=2*pi*f_2
w_3=2*pi*f_3
w_4=2*pi*f_4

w_0=sqrt(w_1*w_2)
f_0=w_0/(2*pi)
W_p=1
W_s=(w_4-w_3)/(w_2-w_1)
bw=w_2-w_1

%% Calculate System's Rank
n=log10(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(log10(W_s))
n=ceil(n)

%% 3dB Frequency
W_0=1/(10^(a_max/10)-1)^(1/(2*n))

%% Butterworth Angles
psi_k=zeros(n,1);
for k=1:n
    psi_k(k)=(pi*(2*k-1)/(2*n)-pi/2);
end

rad2deg(psi_k)

%% Calculate poles
p_k=W_0*(-cos(psi_k)+i*sin(psi_k))'

%% Geffe Algorithm
q_c=w_0/bw
w_0k=zeros(2,ceil(n/2));

for i=1:n % real poles
    if(psi_k(i)==0)
        S_k(i)=-p_k(i);
        Q(i)=q_c/S_k(i);
        psi(i)=acos(1/(2*Q(i)));
        rad2deg(psi(i));
        w_0k(1,i)=w_0;
        poles_i(i)=w_0*(-cos(psi(i))+i*sin(psi(i)))'
    elseif(psi_k(i)<0) % Complex poles
        S_k(i)=-real(p_k(i));
        W_k(i)=-imag(p_k(i));
        C(i)=S_k(i)^2+W_k(i)^2;
        D(i)=2*S_k(i)/q_c;
        E(i)=4+C(i)/q_c^2;
        G(i)=sqrt(E(i)^2-4*D(i)^2);
        Q(i)=1/D(i)*sqrt((E(i)+G(i))/2);
        K(i)=S_k(i)*Q(i)/q_c;
        W(i)=K(i)+sqrt(K(i)^2-1);
        w_0k(1,i)=W(i)*w_0;
        w_0k(2,i)=1/W(i)*w_0;
    end
end
S_k
W_k
C
D
E
G
Q
k
W
w_0k
w0_new=reshape(w_0k.',1,[])';
Q_new=cat(2,Q,Q)';
for i=1:length(w0_new)
    if(w0_new(i)==0)
        Q_new(i)=0;
    end
end
w0_new=nonzeros(w0_new)
Q_new=nonzeros(Q_new)
%% Band Pass Delyiannis-Fried (1st Strategy)

N=cell(n,1);
D=cell(n,1);
T=cell(n,1);

for i=1:n
    R_1old(i)=1;
    C_old(i)=1/(2*Q_new(i))
    R_2old(i)=4*Q_new(i)^2
    
    k_f(i)=w0_new(i)
    k_m(i)=C_old(i)/(k_f(i)*C_avail)
    
    C(i)=C_avail
    R_1(i)=R_1old(i)*k_m(i)
    R_2(i)=R_2old(i)*k_m(i)
    
    N{i}=[2*Q_new(i)*w0_new(i),0];
    D{i}=[1, w0_new(i)/Q_new(i) ,w0_new(i)^2];
    T{i}=tf(N{i},D{i});
    g(i)=abs(freqresp(T{i},w_0));
end

for i=1:n
    plot_transfer_function(T{i},[f_0,f_1,f_2,f_3,f_4])
    name = ['pics/T',num2str(i),'.png'];
    saveas(gcf,name);
end


% for i = 1:n
%     a(i)=1/g(i)
%     Z_2(i)=1/a(i)
%     Z_3(i)=1/(1-a(i))
% end
% gain=10^(LowFreqGain/20)*prod(a)

gain=-10^(LowFreqGain/20)/prod(g)
T_BP=tf(gain,1);
for i=1:n
    T_BP=T{i}*T_BP;
end

plot_transfer_function(T_BP,[f_0,f_1,f_2,f_3,f_4])
saveas(gcf,'pics/T_BP.png');
    

figure('Position', get(0, 'Screensize'));
for i=1:n
    bodemag(T{i});
    hold on;
end
bodemag(T_BP);
grid on;
legend('Unit 1','Unit 2','Unit 3','Unit 4','T_{BP}')
saveas(gcf,'pics/bodeALL.png');

plot_transfer_function(inv(T_BP),[f_0,f_1,f_2,f_3,f_4])
saveas(gcf,'pics/invBP.png');

plot_transfer_function(1/10^(LowFreqGain/20)*T_BP,[f_0,f_1,f_2,f_3,f_4])
saveas(gcf,'pics/T_BP(zero_gain).png');

%% Spectrum Calculation of Input
f1=(w_0-(w_0-w_1)/2)/(2*pi)
f2=(w_0+(w_0+w_1)/3)/(2*pi)
f3=0.4*w_3/(2*pi)
f4=2.5*w_4/(2*pi)
f5=3*w_4/(2*pi)

input= @(t) cos(2*pi*f1*t)+0.8*cos(2*pi*f2*t)+0.8*cos(2*pi*f3*t)+0.6*cos(2*pi*f4*t)+0.5*cos(2*pi*f5*t);
spectrum(T_BP,input)

idealFundFreq = (calcFundamentalFreq([f1;f2;f3;f4;f5]))

function fundFreq = calcFundamentalFreq(F)

F = round(F);
temp = F(end);
for i=length(F)-1:-1:1
    fundFreq = gcd(F(i),temp);
    temp = fundFreq;
end

end

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