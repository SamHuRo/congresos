A = abs(a) ;


lambda0   =   1490e-9;
lambdaRes =   1490e-9 ;

wRes  = 2*pi*c/lambdaRes;
w0    = 2*pi*c/lambda0;

dw  = w0-wRes;


cphia = cumsum( (-wRes*( 1*n2*A.^2/tauR/Seff/n0)*1))*dt;


% cphia =  cphia + cumsum( -wRes*( sigmaN*Np/n0) + dw)*dt;



% cphia(5e5:end) = cphia(5e5);


cphia(1) = cphia(2);

plot(t, cphia)

save cphia cphia