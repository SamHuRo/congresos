% --------------------------------------------------------------------- %
% Self-phase modulation through Kerr effect in a Silicon nanowire dependent
% on the input beam chirp with a Sinc-shaped temporal pulse. Top-Hat
% D-Scan.
% Authors: Samuel Serna - Nicolas Dubreuil
% Date: 18/09/2015
% --------------------------------------------------------------------- %
% clc
% % close all
% clear all
tic,
%% Input parameters

N0 = 6;                             % Pulse duration in number of points
Ntot = 1000;                        % Total number of points for the time frame: EVEN number
lambda = 1.57*10^-4;                % Operating wavelength in cm
c = 3*10^10;                        % speed of light in cm/s
T0 = 1.12*10^-12;                   % Pulse duration in s
phi2max = 3;
n2 = 6*10^-5/100;                       % Nonlinear refractive Index of Silicon in cm^2/GW
beta_TPA = 0.8*0;                     % TPA coefficient in cm/GW
Length = 0.5;                            % Waveguide length in cm
k0 = 2*pi/lambda;                   % wavevector in cm^-1

%% Temporal pulse definition

n = (1:Ntot);
%Pulse0 = sinc((n - Ntot/2)/N0);    % Temporal pulse Careful!!! In mathematica sinc(z)=sin(z)/z. Here sinc(t)=sin(pi*t)/(pi*t) 
Pulse0 = sin(pi.*(n - Ntot/2)/N0)./(pi.*(n - Ntot/2)/N0); 

freq = N0*n/(T0*Ntot);              % in 1/s
dlam = lambda.^2.*freq./c;          % in cm
time = T0*(n-Ntot/2)/N0;

SpectrePulse0 = fftshift(fft(Pulse0));
S0 = (abs(SpectrePulse0)).^2;
count_a = 0;

%% For loop with different chirp values
for PhiNL0 = 0:0.005*pi:1*pi
    count_a = count_a+1;
    I0(count_a) = PhiNL0/(k0*n2*Length);
    TPA_I0_L = beta_TPA*I0(count_a)*Length;
    PhiNL(count_a) = k0*n2/beta_TPA*log(1+TPA_I0_L);
    if beta_TPA == 0
        PhiNL(count_a) = PhiNL0;
    end
    % TPA_I0_L = exp(PhiNL*beta_TPA/(k0*n2))-1;      % beta_TPA*I*L factor uniteless
    % I0(count_a) = TPA_I0_L/(beta_TPA*Length);                     % Initial intensity in GW/cm^2
    depletion = sqrt(I0(count_a)/(1+TPA_I0_L));
    count_b = 0;
    for k = -phi2max:0.05:phi2max
        count_b = count_b+1;
        phi2(count_b) = k;
        L(count_b,:) = exp(1i.*(n-Ntot/2).^2.*k/(2*T0^2.*10^24).*(2*pi*N0/Ntot)^2);
        Pulse(count_b,:) = ifft(SpectrePulse0.*L(count_b,:));
        if depletion == 0
            depletion = 1;
        end
        ChirpedPulse(count_b,:) = Pulse(count_b,:)*depletion.*exp(1i.*PhiNL(count_a).*(abs(Pulse(count_b,:))).^2);
        ChirpedPulsedB(count_b,:) = 10*log10(abs(ChirpedPulse(count_b,:))./max(abs(ChirpedPulse(count_b,:))));
        for j = 1:length(ChirpedPulsedB(count_b,:))
           if ChirpedPulsedB(count_b,j) < -30
               ChirpedPulsedB(count_b,j) = -30;
           end
        end        
        SpectreChirpedPulse(count_b,:) = fft(ChirpedPulse(count_b,:));
        
        S1(count_b,:) = SpectreChirpedPulse(count_b,:).^2;
        NSPI(count_b) = abs(S1(count_b,Ntot/2))/abs(S1(1,Ntot/2));
        S1dB(count_b,:) = 10*log10(abs(S1(count_b,:))./max(abs(S1(count_b,:))));
        for j = 1:length(S1dB(count_b,:))
           if S1dB(count_b,j) < -50
               S1dB(count_b,j) = -50;
           end
        end
        m1(count_b) = sum(dlam.*abs(S1(count_b,:)))/sum(abs(S1(count_b,:)));
        sigma2_1(count_b) = 10^7*(sqrt(sum(dlam.^2.*abs(S1(count_b,:)))/sum(abs(S1(count_b,:))) - m1(count_b).^2));
    end
    sigma2_tot(count_a,:) = sigma2_1;
    delta2sigma2(count_a) = max(sigma2_1)-min(sigma2_1);
    ChirpedPulse_tot(count_a,:,:) = ChirpedPulse;
    ChirpedPulsedB_tot(count_a,:,:) = ChirpedPulsedB;
    SpectreChirpedPulse_tot(count_a,:,:) = SpectreChirpedPulse;
    S_tot(count_a,:,:) = S1;
    SdB_tot(count_a,:,:) = S1dB;
    NSPI_tot(count_a,:) = NSPI;
    PhiNL0_tot(count_a) = PhiNL0;
end
%%
createfigure2sigma(phi2,2.*sigma2_tot')
createfigureSPP(phi2,NSPI_tot')
createfiguredeltasigma(PhiNL./pi,2.*delta2sigma2)
lamver=lambda+dlam-max(dlam)/2;
figure; surf(lamver.*10^7,phi2,S1dB(:,:))
shading interp
figure; surf(time(300:700).*10^12,phi2,ChirpedPulsedB(:,300:700))
shading interp
toc