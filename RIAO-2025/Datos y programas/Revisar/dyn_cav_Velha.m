clear all;
% close all;

n0      =   3.5;              % Linear refractive index
L       =   2*1.6612e-6;      % Cavity length
c       =   3e8;



lambda0 =   1490.00e-9;       % Wavelength
k0      =   2*pi/lambda0;  % Wave vector
tauR    =   n0*L/c;           % Roundtrip time

phi0    =   -.00*pi;           % Phase shift



r       = sqrt(.97);%sqrt(.989);                 % Couping coefficient
k       = sqrt(1-r^2);

a       = .9988;


tau0   = tauR*sqrt(a)/(1-a);
taue    = tauR*sqrt(r)/(1-r);

tau     = 1/(1/taue+1/tau0);

Q       = c*k0*tau/2;

gamma   = 1/tau;

aL      =   -2*log(a)/L;    % Absorption linéaire


T0      =    6e-12;       % Pulse duration 
P0      =    .6;               % Input power
Seff    =    (300e-9)^2;           % Effective mode area

n2      =  6e-18;             % Nonlinear index
beta    = .8e-2/1e9;
sigmaN  = -1.35e-21*1e-6;
sigmaA  = 1.45e-17*1e-4;

hbar    = 6.62606896e-34/(2*pi);

Tmax    = 20*T0;               
M       = 10;                %

dt      = tauR/M;             % Integration step

nt      = 2*Tmax/dt;          % 

t       = (-nt/2:nt/2-1)*dt;  % Time vector
w       = (pi/Tmax)*[(-nt/2:-1) (0:nt/2-1)];

% Init.

Eo      = zeros(1, length(t));
Er      = zeros(1, length(t));
EN      = zeros(1, length(t));
phi     = zeros(1, length(t));
phiE    = zeros(1, length(t));

I       = zeros(M + 2, length(t));
phiNL   = I;
N       = I;
% dn      = Eo;

Chirp_w = unwrap(angle((r - a*exp(1i*w*tauR + 1i*phi0))./(1-r*a*exp(1i*w*tauR + 1i*phi0))));


zero = find(t > -dt & t < dt, 1 );

% atri = t(zero) - T0;
% btri = t(zero);
% ctri = t(zero) + T0;
% 
% tri  = max(min((t-atri)/(btri-atri), (ctri-t)/(ctri-btri)),0);
%rangle = load('angle');

%load('pth');

Ei      = sqrt(P0/Seff)*(exp(-t.^2/(2*T0^2)));%%.*exp(-1i*pth);%
%      Ei_w  = fftshift(fft(Ei)); 
%      Ei_w  = Ei_w.*exp(1*1i.*Chirp_w);
% % 
% %     % Impulsion d'entrée en temporel
%      Ei    = ifft(fftshift(Ei_w));% + sqrt(4*P0)*ones(1,N);

% Loop, all is lossless, in order to reproduce Fig.3 (i)

dz      = dt*L/tauR;

for i = M+1 : nt - M
 
    I(1, i - M) = abs(Er(i-M))^2;
    phiE(i)        = angle(Er(i-M));
 

    %     
     for j = 1 : M+1
%         
        N(j, i - M + j) = N(j, i - M + j - 1) + beta*dt/(2*hbar*c*k0)*(I(j, i - M))^2;
%         
%        N(j, i - M + j) = N(j, i - M + j - 1) + aL*dt/(c*k0)*I(j, i - M);
%         
        I(j + 1, i - M + j) = I(j, i - M + j -  1)*exp(- ( 1*aL + 0*beta*I(j, i - M + j - 1) + 0*sigmaA*N(j, i - M + j - 1) )*dz);
%         
        phiNL(j + 1, i - M + j) = phiNL(j, i - M + j - 1) + dz*k0*( 0*n2*I(j, i - M + j - 1) + 1*sigmaN*N(j, i - M + j - 1) );
%         
     end

    
    EN(i) = exp(-1i*(phi0 + 1*phiNL(j,i)))*sqrt(I(j,i))*exp(1i*phiE(i));

    % Coupler (lossless)
    Eo(i) = r*Ei(i) - 1i*k*EN(i);   
    Er(i) = -1i*k*Ei(i) + r*EN(i);
    
end

phiRes = (sigmaN*N(end-1,:))*k0*L + phi0;

figure(1)

plot(t/tau, abs(Er).^2/P0*Seff, t/tau, abs(Ei).^2/P0*Seff)%, t, trace.trace)

sum(abs(Er).^2)/sum(abs(Ei).^2)

axis([-6 6 0 max(abs(Er).^2)/P0*Seff]);

    aangle = angle(Er);
    

figure(2)

plot(t,  unwrap(aangle)+pi/2, t, unwrap(angle(Ei)))

% axis([-6e-11 6e-11 -pi pi])

% figure(5)
% 
% plot(t, unwrap(angle(Er)));
% 
% 
% 
% axis([-5*T0 5*T0 -2*pi 2*pi]);
% 
% 
% % figure(5)
% % 
% % 
% % dphi = diff(unwrap(angle(Eo)))/(2*pi);
% % 
% %     plot(t(1:end-1), dphi/dt)
% %     
% %     
% %     left = find(t > -5*T0-dt & t < -5*T0+dt);
% %     
% %     right = find(t < 5*T0+dt & t > 5*T0-dt);
% %     
% % axis([-5*T0 5*T0 min(dphi(left(1):right(1))/dt) max(dphi(left(1):right(1))/dt)]);
% 
% figure(2)
% 
% plot(t, abs(Eo).^2);
% 
% axis([-5*T0 5*T0 0 max(abs(Eo).^2)]);
% 
% figure(3)
% 
% plot(abs(Ei).^2, abs(Er).^2)
% 
%  figure(4)
%  
%  A =  P0/Seff*(1 - abs((r - a*exp(1i*w*tauR + 1i*phi0))./(1-r*a*exp(1i*w*tauR + 1i*phi0)))).^2;
%  
%  plot(w, 10*log10(fftshift(abs(fft(fftshift(Er))).^2)), ...); %/max(abs(fft(Eo)).^2)))
%       w, 10*log10(fftshift(abs(fft(fftshift(Ei))).^2)), ...
%       w, 10*log10(A));
%  
%  axis([-2e14 2e14 10*log10(1e-30) max(10*log10(A))]);    
% 