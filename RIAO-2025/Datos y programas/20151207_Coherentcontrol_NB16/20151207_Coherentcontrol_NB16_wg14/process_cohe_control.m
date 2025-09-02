%function process()

nfiles = 313;
sigma2=zeros(nfiles,1);
    for i = 1:nfiles
% 
%         if i < 10 
%             mz = strcat(num2str(0),num2str(0));
%         else
%             mz =num2str(0);
%         end
%         if i>99
%             mz= '';
%         end
%         fn = strcat('W0', mz, num2str(i), '.TXT')
%         importfile(fn);
%         
%         data = evalin('base', 'data');
% 
%         evalin('base', strcat('dat(', num2str(i),', :, :) = data;'));
%         if i <= 9 
%             mz = num2str(0);
%         else
%             mz = '';
%         end
% 
%         fn = strcat('W00', mz, num2str(i), '.txt')
        fn = strcat(num2str(i), '.txt')
        importfile(fn);
        
        data = evalin('base', 'data');

        evalin('base', strcat('dat(', num2str(i),', :, :) = data;'));


    end
    
    Pout = zeros(nfiles, 1);
  
    dat = evalin('base', 'dat');
    Res = .5;
   
    
   for i = 1:nfiles
  
       Span = dat(i,end,1) - dat(i,1,1);
       N    = length(dat(i,:,1));
       spectralin (i,:) = 10.^(dat(i, :, 2)/10);
       Pout(i) = sum(10.^(dat(i, :, 2)/10))*Span/Res/N;
       
                               dlambda = dat(i, :, 1);
                        s = 10.^(dat(i, :, 2)/10);
                        NSPP(i) = s(1,500);
                        m(i) = sum(dlambda.*abs(s))/sum(abs(s));
                        sigma2(i) = sqrt(sum(dlambda.^2.*abs(s))/sum(abs(s)) - m(i).^2);
       
       
       figure(1); hold on
       
       for j = 1:length(dat(i, : ,2))
           if dat(i, j ,2) < -80
               dat(i, j ,2) = -80;
           end
       end
       
       %plot3(dat(i, : ,1), i*ones(1,length(dat(i,:,1))), -(dat(i, : ,2))/max((dat(i, : ,2))))
   end
   
%                   view(0,70)
%    print('-dpng', 'spectres.png');
%    
%    xlabel('$$\lambda$$ [nm]', 'interpreter', 'latex');
%    ylabel('Spectre nb', 'interpreter', 'latex');
%    zlabel('Power [dB]', 'interpreter', 'latex');
%   
%   figure(i+1)
%    
%    assignin('base', 'Pout', Pout);
%    
%    %Pin = load('power.txt');
%     Pin=[0.150 1 2.5 4 6];  
%       
%       assignin('base', 'Pin', Pin);
%       
%       %plot(Pin.^2, Pin.^2./Pout.^2, 'o')
%       %plot(Pin, Pin./Pout, 'o')
%       plot(Pin,Pout, 'o')
%       xlabel('Pin [mW]', 'interpreter', 'latex');
%       ylabel('Pout [mW]', 'interpreter', 'latex');
%       
%           print('-dpng', 'PinPout.png');
%       figure;
%       
%             plot(Pin,Pin'./Pout, 'o')
%       xlabel('Pin [mW]', 'interpreter', 'latex');
%       ylabel('Pin/Pout', 'interpreter', 'latex'); 
%       
%          print('-dpng', 'PinPout_fPin.png');
%       
%                            hf = figure;
%                         hold on;
%                         plot(Pin, sqrt(sigma2), 'bo');
%                         
%                                      print('-dpng', 'sigma2.png');  
%                                                                         
%                         set(hf, 'Position', [900 400 800 700])
%                         ylabel('$$\sigma\,$$[{nm}\,]', 'interpreter', 'latex')
%                         xlabel('Input power [mW]', 'interpreter', 'latex')
%end

%%
Position2=0:25;
Phi2ps=zeros(1,26);
% Phi2ps(1,k+1)=0.228.*(Position(1,k+1))-2.83;
for k=0:25
    Phi2ps(1,k+1)=0.2251.*(Position2(1,k+1))-2.75+0.1;
end

plot(Phi2ps,2*sigma2(287:312,1),'LineWidth',2);
hold on
plot(Phi2ps,2*sigma2(261:286,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(235:260,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(209:234,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(183:208,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(157:182,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(131:156,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(105:130,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(79:104,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(53:78,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(27:52,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(1:26,1),'LineWidth',2);
hold on
% axis ([-3.8 4.2 3 12])
% xlabel('Stage position [mm]', 'interpreter', 'latex')
xlim([-3 3])
box on
grid on
set(gca,'FontSize',12,'linewidth',2)
xlabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',14,'interpreter', 'latex')
ylabel('$2\sigma$ [nm]','fontsize',14, 'interpreter', 'latex')
% legend('P_{in}=10 mW','P_{in}=7 mW','P_{in}=2.5 mW','P_{in}=1 mW','P_{in}=0.1 mW')
%% Baricenter
figure
plot(Phi2ps,m(287:312),'LineWidth',2);
hold on
plot(Phi2ps,m(261:286),'LineWidth',2);
plot(Phi2ps,m(235:260),'LineWidth',2);
plot(Phi2ps,m(209:234),'LineWidth',2);
plot(Phi2ps,m(183:208),'LineWidth',2);
plot(Phi2ps,m(157:182),'LineWidth',2);
plot(Phi2ps,m(131:156),'LineWidth',2);
plot(Phi2ps,m(105:130),'LineWidth',2);
plot(Phi2ps,m(79:104),'LineWidth',2);
plot(Phi2ps,m(53:78),'LineWidth',2);
plot(Phi2ps,m(27:52),'LineWidth',2);
plot(Phi2ps,m(1:26),'LineWidth',2);
hold on
% axis ([-3.8 4.2 3 12])
% xlabel('Stage position [mm]', 'interpreter', 'latex')
xlim([-3 3])
box on
grid on
set(gca,'FontSize',12,'linewidth',2)
xlabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',14,'interpreter', 'latex')
ylabel('Barycenter [nm]','fontsize',14, 'interpreter', 'latex')
% legend('P_{in}=10 mW','P_{in}=7 mW','P_{in}=2.5 mW','P_{in}=1 mW','P_{in}=0.1 mW')
%% 0.1mW

        figure;
       for ms=1:26
       plot3(dat(ms, : ,1), Phi2ps(1,ms)*ones(1,length(dat(i,:,1))), (10^6.*spectralin(ms, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^6.*spectralin(313, : ,1)))
       view(0,70)
             
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 1mW  
       figure;
       for mp=27:52           
       plot3(dat(mp, : ,1), Phi2ps(1,(mp-26))*ones(1,length(dat(mp,:,1))), (10^6.*spectralin(mp, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
       view(0,70)
       
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 2mW  
   
       figure;
       for mt=53:78           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-52))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
       view(0,70)
       
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 3mW  
   
       figure;
       for mt=79:104           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-78))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)       

   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 4mW  
   
       figure;
       for mt=105:130           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-104))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 5mW  
   
       figure;
       for mt=131:156           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-130))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 6mW  
   
       figure;
       for mt=157:182           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-156))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 7mW  
   
       figure;
       for mt=183:208           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-182))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 8mW  
   
       figure;
       for mt=209:234           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-208))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 9mW  
   
       figure;
       for mt=235:260           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-234))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');
%% 12mW  
   
       figure;
       for mt=261:286           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-260))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');

%% 15mW  

       figure;
       for mt=287:312           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-286))*ones(1,length(dat(mt,:,1))), (10^6.*spectralin(mt, : ,1)))
       hold on
       end
       hold on
       plot3(dat(313, : ,1), (Phi2ps(1,26)+0.1)*ones(1,length(dat(i,:,1))), (10^7.*spectralin(313, : ,1)))
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [nW]','fontsize',13,'interpreter', 'latex');
    xlim([1580 1590])
    ylim([-3 3.2])
%    zlim([-80 -20])
%    print('-depsc2', 'Pin0p1mWSpectra.eps');

   
%%


figure
PoutnW=Pout.*1000000;
plot(Phi2ps,PoutnW(287:312,1),'LineWidth',2);
hold on
plot(Phi2ps,PoutnW(261:286,1),'LineWidth',2);
plot(Phi2ps,PoutnW(235:260,1),'LineWidth',2);
plot(Phi2ps,PoutnW(209:234,1),'LineWidth',2);
plot(Phi2ps,PoutnW(183:208,1),'LineWidth',2);
plot(Phi2ps,PoutnW(157:182,1),'LineWidth',2);
plot(Phi2ps,PoutnW(131:156,1),'LineWidth',2);
plot(Phi2ps,PoutnW(105:130,1),'LineWidth',2);
plot(Phi2ps,PoutnW(79:104,1),'LineWidth',2);
plot(Phi2ps,PoutnW(53:78,1),'LineWidth',2);
plot(Phi2ps,PoutnW(27:52,1),'LineWidth',2);
plot(Phi2ps,PoutnW(1:26,1),'LineWidth',2);
xlim([-3 3])
box on
grid on
set(gca,'FontSize',12,'linewidth',2)
xlabel('$\phi^{(2)} (\rm ps^{2})$','fontsize',14,'interpreter', 'latex')
ylabel('Pout (nW)','fontsize',14, 'interpreter', 'latex')
%legend('P_{in}=10 mW','P_{in}=7 mW','P_{in}=2.5 mW','P_{in}=1 mW','P_{in}=0.1 mW','Location','NorthEastOutside')


%%
[X,Y] = meshgrid(Phi2ps,data(:,1));
spectralin0p1mW=spectralin(1:26,:);
surfpower(Y,X,spectralin0p1mW'.*10^6)

%%    
spectralin1mW=spectralin(27:52,:);
surfpower(Y,X,spectralin1mW'.*10^6)
%%    
spectralin2mW=spectralin(53:78,:);
surfpower(Y,X,spectralin2mW'.*10^6)
%%    
spectralin3mW=spectralin(79:104,:);
surfpower(Y,X,spectralin3mW'.*10^6)
%%    
spectralin4mW=spectralin(105:130,:);
surfpower(Y,X,spectralin4mW'.*10^6)
%%    
spectralin5mW=spectralin(131:156,:);
surfpower(Y,X,spectralin5mW'.*10^6)
%%    
spectralin6mW=spectralin(157:182,:);
surfpower(Y,X,spectralin6mW'.*10^6)
%%    
spectralin7mW=spectralin(183:208,:);
surfpower(Y,X,spectralin7mW'.*10^6)
%%    
spectralin8mW=spectralin(209:234,:);
surfpower(Y,X,spectralin8mW'.*10^6)
%%    
spectralin9mW=spectralin(235:260,:);
surfpower(Y,X,spectralin9mW'.*10^6)
%%    
spectralin12mW=spectralin(261:286,:);
surfpower(Y,X,spectralin12mW'.*10^6)
%%    
spectralin15mW=spectralin(287:312,:);
surfpower(Y,X,spectralin15mW'.*10^6)
% figure
% PoutuW=Pout.*1000;
% Pin2=[0.1 1 2.5 7 10];
% [X,Y] = meshgrid(Phi2ps,Pin2);
% PoutTot(1,:)=ones(size(Phi2ps)).*PoutuW(13,1);
% PoutTot(2,:)=PoutuW(27:55,1);
% PoutTot(3,:)=PoutuW(56:84,1);
% PoutTot(4,:)=PoutuW(85:113,1);
% PoutTot(5,:)=PoutuW(114:142,1);
% surf(X,Y,PoutTot)
% shading interp
% xlabel('$\phi^{(2)} (\rm ps^{2})$','fontsize',14,'interpreter', 'latex')
% ylabel('Pin (mW)','fontsize',14, 'interpreter', 'latex')
% zlabel('P_{out}(\muW)')
% view(-132.5,35)
% colorbar
%%
% figure
% PinPoutTot(1,:)=Pin2(1)./PoutTot(1,:).*1000;
% PinPoutTot(2,:)=Pin2(2)./PoutTot(2,:).*1000;
% PinPoutTot(3,:)=Pin2(3)./PoutTot(3,:).*1000;
% PinPoutTot(4,:)=Pin2(4)./PoutTot(4,:).*1000;
% PinPoutTot(5,:)=Pin2(5)./PoutTot(5,:).*1000;
% surf(X,Y,PinPoutTot)
% shading interp
% xlabel('$\phi^{(2)} (\rm ps^{2})$','fontsize',14,'interpreter', 'latex')
% ylabel('Pin (mW)','fontsize',14, 'interpreter', 'latex')
% zlabel('P_{in}/P_{out}')
% view(-132.5,35)
% colorbar
