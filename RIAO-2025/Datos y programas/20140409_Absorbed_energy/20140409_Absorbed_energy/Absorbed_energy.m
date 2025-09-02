%%
importdata('data.mat');
n0=3.48;
c=3e8;
lambda0=1578.8e-9;
VFCR=1.7*(lambda0/n0)^3;
hbar=1.054e-34;
omega0=2*pi*c/lambda0;
neff=2.71;
betaTPA=8.4e-12;
SigmaR=1.34e-27;

deltat=5.8e-12;

deltalambda=1.7e-9;
Energy0=n0./c.*VFCR.*sqrt(2*hbar*omega0*neff/(betaTPA*SigmaR*lambda0)*deltalambda./deltat);
Energy=n0./c.*VFCR.*sqrt(2*hbar*omega0*neff/(betaTPA*SigmaR*lambda0)*yy2.*1e-9./deltat);
figure(20)
   plot(phi2,Energy.*1e12,'LineWidth',2)
   grid on
   %axis ([-3 3.5 1.2 3.5])
   xlabel('$\varphi^{(2)} (\rm ps^{2})$', 'interpreter', 'latex','FontSize',20);
   ylabel('$ \rm Energy [\rm pJ]$', 'interpreter', 'latex','FontSize',20);
   %plotTickLatex2D('FontSize',20)
   axis([-3 3.5 0.165 0.28])