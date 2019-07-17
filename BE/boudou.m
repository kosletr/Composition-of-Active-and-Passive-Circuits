format long
%%
%stoixeia aem
%%
a1 = 8;
a2 = 8;
a3 = 5;
a4 = 1;
%%
%prodiagrafes
%%
f0 = 2400

f1 = 1725 + 25*a4

f2 = (f0^2) / f1

D = (1 / 2.2) * ((f0^2 - f1^2) / f1) 

f3 = (-D + sqrt(D^2 + 4*f0^2)) / 2

f4 = (f0^2) / f3

amin = 28 + a3 * (5/9)

amax = 0.5 + (a4 / 18)
%%
%loipa
%%
omega1 = 2 * pi * f1 %ù1

omega2 = 2 * pi * f2 %ù2

omega3 = 2 * pi * f3 %ù3

omega4 = 2 * pi * f4 %ù4

omega0 = sqrt(omega1 * omega2) %ù0

bw = omega2 - omega1 %bandwidth

OMEGAP = 1 %Ùp

OMEGAS =  ( omega2 - omega1 ) / (omega4 - omega3) %Ùs

n = acosh( sqrt( (((10^(amin/10)) -1) / ( (10^(amax/10)) -1 ))) ) / acosh(OMEGAS)
%epilegoyme ton amesos epomeno akeraio: dinadi n = 4

n= ceil(n)

e = 1 /( sqrt ( ((10^(amin/10)) -1) ) ) %å

a = (1/n) * asinh(1/e) % á

%%
% gonies butterworth gia n = 4
angle0 = pi / 8; %ø1,2 = 22.5

angle1 = 3 * pi / 8;  %ø3,4 = 67.5
%% etsi epeidi mporo
koa = cosh(a) 

sia = sinh(a)
%%
%upologismos polwn: pragmatiko kai fantasiko meros
%% sinolo polwn: 4
%%
real1 = sia * cos(angle0)

imagine1 = koa* sin(angle0)

real2 = sia * cos(angle1)

imagine2 = koa* sin(angle1)
%%
OMEGA012 = sqrt ( real1^2 + imagine1^2 ) %Ù01,2

OMEGA034 = sqrt ( real2^2 + imagine2^2 ) %Ù03,4

Q12 = OMEGA012 / (2 * real1) 

Q34 = OMEGA034 / (2 * real2)

%%
%poloi poy prokyptoyn apo tin antistrofi ton polon chebyshev
%%

tOMEGA012 = 1 / OMEGA012 

tOMEGA034 = 1 / OMEGA034
%%
%klimakopoioume ta metra ton polon oste na metaferthoume sto pedio ton
%sixnotiton
%%
tkOMEGA012 = OMEGAS * tOMEGA012
tkOMEGA034 = OMEGAS * tOMEGA034

%ta midenika tis apokrisis ICH tha einai:
OMEGAZ1 = sec (pi / (2*n) )
OMEGAZ2 = sec (3*pi / (2*n)) 
%%
%klimakopoioyme ta midenika tis ICH
%%
tOMEGAZ1 = OMEGAS * OMEGAZ1
tOMEGAZ2 = OMEGAS * OMEGAZ2
%%
%antistrefoyme toys poloys tis ICH
%%
hOMEGA012 = 1 / tkOMEGA012
hOMEGA034 = 1 / tkOMEGA034
%%
%antistrefoume kai ta midenika tis ICH
%%
hOMEGAZ1 = 1 / tOMEGAZ1
hOMEGAZ2 = 1 / tOMEGAZ2
%%
%oi poloi tis anodiabatis sinartisis einai
hS12 =  (-hOMEGA012) / (2*Q12)
hS34 = (-hOMEGA034) / (2*Q34)

hOMEGA12 = sqrt( (hOMEGA012)^2 - (hS12)^2 )

hOMEGA34 = sqrt( (hOMEGA034)^2 - (hS34)^2 )
%%
qc = omega0 / bw 
%% akolouthei metasximatismos migadikoy polou
real12 = -hS12
imagine12 = -hOMEGA12
C12 = real12^2 + imagine12^2
D12 = (2*real12) / qc
E12 = 4 + C12 / (qc*qc)
G12 = sqrt(E12^2 - 4*D12^2)
Q1_2 = (1/D12) * sqrt((1/2) *(E12+G12))
k12 = (real12* Q1_2) / qc
W12 = k12 + sqrt(k12^2-1)
omega01 = W12 * omega0
omega02 =(1/ W12) * omega0

%% akolouthei metasximatismos migadikoy polou
real34 = -hS34
imagine34 = -hOMEGA34
C34 = real34^2 + imagine34^2
D34 = (2*real34) / qc
E34 = 4 + C34 / (qc*qc)
G34 = sqrt(E34^2 - 4*D34^2)
Q3_4 = (1/D34) * sqrt((1/2) *(E34+G34))
k34 = (real34* Q3_4) / qc
W34 = k34 + sqrt(k34^2-1)
omega03 = W34 * omega0
omega04 =(1/ W34) * omega0

%% epomeni kinisi: metasximatismos midenikou
K = 2+ hOMEGAZ1^2 /  qc^2
x = (K + sqrt(K^2-4)) / 2
omegaz1 = omega0 * sqrt(x)
omegaz2 = omega0 * (1 / sqrt(x))
%% methepomeni kinisi: metasximatismos midenikou
K = 2+ hOMEGAZ2^2 /  qc^2
x = (K + sqrt(K^2-4)) / 2
omegaz3 = omega0 * sqrt(x)
omegaz4 = omega0 * (1 / sqrt(x))

t1 = tf([1 0 omegaz1^2],[1 (omega01/Q1_2) omega01^2])

t2 = tf([1 0 omegaz2^2],[1 (omega02/Q1_2) omega02^2])

t3 = tf([1 0 omegaz3^2],[1 (omega03/Q3_4) omega03^2])

t4 = tf([1 0 omegaz4^2],[1 (omega04/Q3_4) omega04^2])

t= t1*t2*t3*t4
plot_transfer_function( t, [f0 f1 f2 f3 f4]);

 

