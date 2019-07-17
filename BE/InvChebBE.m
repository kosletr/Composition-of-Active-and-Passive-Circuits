% Composition of Active and Passive Circuits 2019
% Inverse Chebyshev Band Elimination Filter
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

f_0=2400
f_1=1725+25*a4
f_2=f_0^2/f_1
D=(1/2.2)*(f_0^2-f_1^2)/f_1;
f_3=(-D+sqrt(D^2+4*f_0^2))/2
f_4=f_0^2/f_3
a_min=28+a3*5/9
a_max=0.5+a4/18
C_avail=0.01e-6
LowFreqGain=10


% f_1=1000/(2*pi)
% f_2=4000/(2*pi)
% f_3=1600/(2*pi)
% f_4=2500/(2*pi)
% a_min=25
% a_max=1
% C_avail=0.1e-6
% LowFreqGain=0

% f_1=3998.65/(2*pi)
% f_2=39451.63/(2*pi)
% f_3=8478/(2*pi)
% f_4=18607.38/(2*pi)
% a_min=23.88
% a_max=0.88
% C_avail=0.01e-6
% LowFreqGain=5

%% Calculate BW
w_1=2*pi*f_1
w_2=2*pi*f_2
w_3=2*pi*f_3
w_4=2*pi*f_4

w_0=sqrt(w_1*w_2)
f_0=w_0/(2*pi)
W_p=1
W_s=(w_2-w_1)/(w_4-w_3)
bw=w_2-w_1

%% Calculate System's Rank
n=acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(W_s))
n=ceil(n)
epsilon=(10^(a_min/10)-1)^(-1/2)
a=1/n*asinh(1/epsilon)

%% Butterworth Angles
psi_k=zeros(n,1);
for k=1:n
    psi_k(k)=(pi/2-pi*(2*k-1)/(2*n));
end
psi_k=circshift(psi_k,[2 2]);
rad2deg(psi_k)

%% Calculate Poles
Sigma_k=(sinh(a)*cos(psi_k))';
W_k=(cosh(a)*sin(psi_k))';
Wo_k=(sqrt(Sigma_k.^2+W_k.^2))
Q_k=(Wo_k./(2*Sigma_k))

p_k=(-Sigma_k+i*W_k)'
tilde_Wo_k=W_s./Wo_k

%% Calculate Zeros
for k=1:floor(n/2)
    W_z(k)=sec((2*k-1)*pi/(2*n));
end
if(ceil(n/2)~=n/2)
    W_z(end+1)=Inf;
end
W_z
tilde_W_z=(W_s.*W_z);
tilde_W_z

%% Calculate Inverse Poles and Zeros
invP_k=(1./tilde_Wo_k)
hat_W_z=1./tilde_W_z

for i=1:n
    if(psi_k(i)==0)
        hat_Sigma_k(i)=-invP_k(i);
    else
        hat_Sigma_k(i)=-invP_k(i)/(2*Q_k(i));
        hat_W_k(i)=sqrt(invP_k(i)^2-hat_Sigma_k(i)^2);
    end
end
hat_Sigma_k
hat_W_k

%% Geffe Algorithm
q_c=w_0/bw
for i=1:floor(n/2)
    K(i)=2+hat_W_z(i)^2/q_c^2;
end
K_fl=fliplr(K);
if(ceil(n/2)~=n/2)
    i=ceil(n/2);
    K=[K,2+hat_W_z(i)^2/q_c^2,K_fl]
else
    K=[K,K_fl]
end

for i=1:n
    if imag(p_k(i))==0
        Q_temp(i)=-q_c/hat_Sigma_k(i);
        psi(i)=acos(1/(2*Q_temp(i)));
        rad2deg(psi(i));
        w_0k(i)=w_0
    else
        C(i)=hat_Sigma_k(i)^2+hat_W_k(i)^2;
        D(i)=-2*hat_Sigma_k(i)/q_c;
        E(i)=4+C(i)/q_c^2;
        G(i)=sqrt(E(i)^2-4*D(i)^2);
        Q_temp(i)=1/D(i)*sqrt((E(i)+G(i))/2);
        k(i)=-hat_Sigma_k(i)*Q_temp(i)/q_c;
        W(i)=k(i)+sqrt(k(i)^2-1);
        x(i)=(K(i)+sqrt(K(i)^2-4))/2;
        w_0k(2*i-1)=W(i)*w_0;
        w_0k(2*i)=1/W(i)*w_0;
        w_z(2*i-1)=w_0*sqrt(x(i));
        w_z(2*i)=w_0/sqrt(x(i));
    end
end
hat_Sigma_k
hat_W_k
C
D
E
G
Q_temp
k
W
K
x
%w_0k=reshape(w_0k.',1,[])
w_0k
w_z
Q=zeros(1,n);
for i=1:floor(n/2)
    Q(2*i-1)=Q_temp(i);
    Q(2*i)=Q_temp(i);
end
if(ceil(n/2)~=n/2)
    Q(ceil(n/2))=Q_temp(ceil(n/2));
end
Q
%% High Pass - Low Pass - Notch

Num=cell(n,1);
Den=cell(n,1);
T=cell(n,1);


for i=1:n
    if(w_z(i)/w_0k(i) <= 1) %HPN & Notch
        %w_z(i)
        %w_0k(i)
        k1(i)=(w_0k(i)/w_z(i))^2-1
        R_1old(i)=1;
        R_3old(i)=1;
        R_2old(i)=(2+k1(i))^2*Q(i)^2
        R_4old(i)=(2+k1(i))*Q(i)^2
        C_old(i)=1/((2+k1(i))*Q(i))
        k2(i)=(2+k1(i))*Q(i)^2/((2+k1(i))*Q(i)^2+1)
        k(i)=k2(i)*(w_0k(i)/w_z(i))^2
        k_f(i)=w_0k(i)
        k_m(i)=C_old(i)/(k_f(i)*C_avail)
        
        C(i)=C_avail
        C1(i)=k1(i)*C_avail;
        R_1(i)=R_1old(i)*k_m(i)
        R_2(i)=R_2old(i)*k_m(i)
        R_3(i)=R_3old(i)*k_m(i)
        R_4(i)=R_4old(i)*k_m(i)
        
        Num{i}=k(i)*[1,(k2(i)-1)/(k2(i)*(k1(i)+1)*R_1(i)*C(i))+(k1(i)+2)/((k1(i)+1)*R_2(i)*C(i)),1/(R_1(i)*R_2(i)*C(i)^2*(k1(i)+1))];
        Den{i}=[1,(k1(i)+2)/(R_2(i)*C(i)),1/(R_1(i)*R_2(i)*C(i)^2)];
        T{i}=tf(Num{i},Den{i});
        g(i)=abs(freqresp(T{i},0));
    else % LPN
        %w_z(i)
        %w_0k(i)
        omega_z(i)=w_z(i)/w_0k(i)
        R_1old(i)=1;
        R_4old(i)=1;
        C_old(i)=1/(2*Q(i))
        R_2old(i)=4*Q(i)^2
        R_5old(i)=4*Q(i)^2/(omega_z(i)^2-1)
        R_3old(i)=omega_z(i)^2/(2*Q(i)^2)
        kH(i)=1/(1+omega_z(i)^2/(2*Q(i)^2))
        kL(i)=kH(i)*omega_z(i)^2
        k_f(i)=w_0k(i)
        k_m(i)=C_old(i)/(k_f(i)*C_avail)
        
        C(i)=C_avail
        R_1(i)=R_1old(i)*k_m(i)
        R_2(i)=R_2old(i)*k_m(i)
        R_3(i)=R_3old(i)*k_m(i)
        R_4(i)=R_4old(i)*k_m(i)
        R_5(i)=R_5old(i)*k_m(i)
        kH(i)=1/(1+omega_z(i)^2/(2*Q(i)^2))
        kL(i)=kH(i)*omega_z(i)^2
        
        Num{i}=k(i)*[1,((k(i)-1)/k(i)*1/(R_1(i)*C(i))+2/(R_2(i)*C(i))+2/(R_5(i)*C(i))),(1/(R_1(i)*R_5(i)*C(i)^2)+1/(R_1(i)*R_2(i)*C(i)^2))];
        Den{i}=[1,2/(R_2(i)*C(i)),1/(R_1(i)*R_2(i)*C(i)^2)];
        T{i}=tf(Num{i},Den{i});
        g(i)=abs(freqresp(T{i},0));
    end
end

for i=1:n
    plot_transfer_function(T{i},[f_0,f_1,f_2,f_3,f_4])
    name = ['pics/T',num2str(i),'.png'];
    saveas(gcf,name);
end

K_g=-10^(LowFreqGain/20)/prod(g)
r1=1e4
r2=-K_g*r1
T_BE=tf(K_g,1);
for i=1:n
    T_BE=T{i}*T_BE;
end

plot_transfer_function(T_BE,[f_0,f_1,f_2,f_3,f_4])
saveas(gcf,'pics/T_BE.png');

figure;
figure('Position', get(0, 'Screensize'));
for i=1:n
    bodemag(T{i});
    hold on;
end
bodemag(T_BE);
grid on;
legend('Unit 1','Unit 2','Unit 3','Unit 4','T_{BE}')
saveas(gcf,'pics/bodeALL.png');

plot_transfer_function(inv(T_BE),[f_0,f_1,f_2,f_3,f_4])
saveas(gcf,'pics/invBP.png');

plot_transfer_function(1/10^(LowFreqGain/20)*T_BE,[f_0,f_1,f_2,f_3,f_4])
saveas(gcf,'pics/checkSpecs.png');

Max_Spec = -a_max+LowFreqGain
Min_Spec = -a_min+LowFreqGain

% testTF=tf(10^(Gain/20),1);
% for i=1:n
%     Nu{i}=[1,0,w_z(i)^2];
%     De{i}=[1,w_0k(i)/Q(i),w_0k(i)^2];
%     testTF=testTF*tf(Nu{i},De{i});
% end
% plot_transfer_function(testTF,[f_0,f_1,f_2,f_3,f_4])
% plot_transfer_function(inv(testTF),[f_0,f_1,f_2,f_3,f_4])

%% Spectrum Calculation of Input
f1=w_0-(w_0-w_3)/2
f2=w_0+(w_0+w_3)/2
f3=0.5*w_1
f4=2.4*w_2
f5=3.5*w_2

input= @(t) 0.8*cos(f1*t)+cos(f2*t)+cos(f3*t)+0.8*cos(f4*t)+0.4*cos(f5*t);
spectrum(T_BE,input)

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