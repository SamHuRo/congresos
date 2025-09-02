%function process()

nfiles = 183;
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
Phi2ps2=zeros(1,26);
% Phi2ps(1,k+1)=0.228.*(Position(1,k+1))-2.83;
for k=0:25
    Phi2ps2(1,k+1)=0.2251.*(Position2(1,k+1))-2.75+0.1;
end
Position=0:28;
Phi2ps=zeros(1,29);
for k=0:9
    Phi2ps(1,k+1)=0.2251.*(Position(1,k+1))-2.75+0.1;
end
for k=10:15
    Phi2ps(1,k+1)=0.2251.*(Position(1,k+1)/2+5)-2.75+0.1;
end
for k=16:28
    Phi2ps(1,k+1)=0.2251.*(Position(1,k+1)-3)-2.75+0.1;
end
plot(Phi2ps,2*sigma2(114:142,1),'LineWidth',2);
hold on
plot(Phi2ps,2*sigma2(85:113,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(56:84,1),'LineWidth',2);
plot(Phi2ps,2*sigma2(27:55,1),'LineWidth',2);
plot(Phi2ps2,2*sigma2(1:26,1),'LineWidth',2);
hold on
% axis ([-3.8 4.2 3 12])
% xlabel('Stage position [mm]', 'interpreter', 'latex')
xlim([-3 3])
box on
grid on
set(gca,'FontSize',12,'linewidth',2)
xlabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',14,'interpreter', 'latex')
ylabel('$2\sigma$ [nm]','fontsize',14, 'interpreter', 'latex')
legend('P_{in}=10 mW','P_{in}=7 mW','P_{in}=2.5 mW','P_{in}=1 mW','P_{in}=0.1 mW')

% Position=0:38;
% Phi2ps=zeros(1,39);
% 
% for k=0:5
%     Phi2ps(1,k+1)=0.228.*(Position(1,k+1))-2.63;
% end
% for k=6:30
%     Phi2ps(1,k+1)=0.228.*(Position(1,k+1)/2+3)-2.63;
% end
% for k=31:38
%     Phi2ps(1,k+1)=0.228.*(Position(1,k+1)-12)-2.63;
% end
% plot(Phi2ps,sigma2(1:39,1),'LineWidth',2);
% plot(Phi2ps,sigma2(40:78,1),'Color','m','LineWidth',2);
% plot(Phi2ps,sigma2(79:117,1),'Color','g','LineWidth',2);
% plot(Phi2ps,sigma2(118:156,1),'Color','k','LineWidth',2);

% sigma2_9mW=zeros(26,1);
% for k=5:5:1305
%     sigma2_9mW(k/5,1)=sigma2(k,1);
% end
% plot(Position,flipud(sigma2_9mW),'LineWidth',2);%,'o');
% hold on
% 
% sigma2_4mW=zeros(26,1);
% for k=4:5:129
%     sigma2_4mW((k+1)/5,1)=sigma2(k,1);
% end
% plot(Position, flipud(sigma2_4mW),'Color','m','LineWidth',2);
% 
% sigma2_2500uW=zeros(26,1);
% for k=3:5:128
%     sigma2_2500uW((k+2)/5,1)=sigma2(k,1);
% end
% plot(Position, flipud(sigma2_2500uW),'Color','g','LineWidth',2);
% 
% sigma2_1mW=zeros(26,1);
% for k=2:5:127
%     sigma2_1mW((k+3)/5,1)=sigma2(k,1);
% end
% plot(Position,flipud(sigma2_1mW),'Color','k','LineWidth',2);
% 
% sigma2_100uW=zeros(26,1);
% for k=1:5:126
%     sigma2_100uW((k+4)/5,1)=sigma2(k,1);
% end
% plot(Position,flipud(sigma2_100uW),'Color','r','LineWidth',2);
%axis ([0 25 3 12])
%plot(Phi2ps,sigma2(1:39,1),'LineWidth',2);%,'o');
%%

        figure;
       for ms=1:26
       plot3(dat(ms, : ,1), Phi2ps2(1,ms)*ones(1,length(dat(i,:,1))), (dat(ms, : ,2)))
       hold on
       end
       view(0,70)
             
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [dBm]','fontsize',13,'interpreter', 'latex');
   xlim([1553 1593])
   ylim([-3 3])
   zlim([-80 -20])
   print('-depsc2', 'Pin0p1mWSpectra.eps');
   
       figure;
       for mp=27:55           
       plot3(dat(mp, : ,1), Phi2ps(1,(mp-26))*ones(1,length(dat(mp,:,1))), (dat(mp, : ,2)))
       hold on
       end
       view(0,70)
       
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [dBm]','fontsize',13,'interpreter', 'latex');
   xlim([1553 1593])
   ylim([-3 3])
   zlim([-80 -20])
   print('-depsc2', 'Pin1mWSpectra.eps');
   
       figure;
       for mt=56:84           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-55))*ones(1,length(dat(mt,:,1))), (dat(mt, : ,2)))
       hold on
       end
       view(0,70)
       
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [dBm]','fontsize',13,'interpreter', 'latex');
   xlim([1553 1593])
   ylim([-3 3])
   zlim([-80 -20])
   print('-depsc2', 'Pin2p5mWSpectra.eps');
   
       figure;
       for mt=85:113           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-84))*ones(1,length(dat(mt,:,1))), (dat(mt, : ,2)))
       hold on
       end
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [dBm]','fontsize',13,'interpreter', 'latex');
   xlim([1553 1593])
   ylim([-3 3])
   zlim([-80 -20])
   print('-depsc2', 'Pin7mWSpectra.eps');

       figure;
       for mt=114:142           
       plot3(dat(mt, : ,1), Phi2ps(1,(mt-113))*ones(1,length(dat(mt,:,1))), (dat(mt, : ,2)))
       hold on
       end
       view(0,70)
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [dBm]','fontsize',13,'interpreter', 'latex');
   xlim([1553 1593])
   ylim([-3 3])
   zlim([-80 -20])
   print('-depsc2', 'Pin10mWSpectra.eps');
   
 %%  
       figure;
       for ms=1:26
       plot3(dat(ms, : ,1), Phi2ps2(1,ms)*ones(1,length(dat(i,:,1))), (10^3.*spectralin(ms, : ,1)))
       hold on
       end
       view(0,70)
             
   xlabel('$\lambda [\rm nm]$','fontsize',13,'interpreter', 'latex');
   ylabel('$\phi^{(2)} [\rm ps^{2}]$','fontsize',13,'interpreter', 'latex');
   zlabel('Power [$\mu$W]','fontsize',13,'interpreter', 'latex');
   xlim([1553 1593])
   ylim([-3 3])
   print('-depsc2', 'Pin0p1mWLinearSpectra.eps');
%%
figure
PoutuW=Pout.*1000;
plot(Phi2ps,PoutuW(114:142,1),'LineWidth',2);
hold on
plot(Phi2ps,PoutuW(85:113,1),'LineWidth',2);
plot(Phi2ps,PoutuW(56:84,1),'LineWidth',2);
plot(Phi2ps,PoutuW(27:55,1),'LineWidth',2);
plot(Phi2ps2,PoutuW(1:26,1),'LineWidth',2);
xlim([-3 3])
box on
grid on
set(gca,'FontSize',12,'linewidth',2)
xlabel('$\phi^{(2)} (\rm ps^{2})$','fontsize',14,'interpreter', 'latex')
ylabel('Pout ($\mu$W)','fontsize',14, 'interpreter', 'latex')
legend('P_{in}=10 mW','P_{in}=7 mW','P_{in}=2.5 mW','P_{in}=1 mW','P_{in}=0.1 mW','Location','NorthEastOutside')


%%
figure
PoutuW=Pout.*1000;
Pin2=[0.1 1 2.5 7 10];
[X,Y] = meshgrid(Phi2ps,Pin2);
PoutTot(1,:)=ones(size(Phi2ps)).*PoutuW(13,1);
PoutTot(2,:)=PoutuW(27:55,1);
PoutTot(3,:)=PoutuW(56:84,1);
PoutTot(4,:)=PoutuW(85:113,1);
PoutTot(5,:)=PoutuW(114:142,1);
surf(X,Y,PoutTot)
shading interp
xlabel('$\phi^{(2)} (\rm ps^{2})$','fontsize',14,'interpreter', 'latex')
ylabel('Pin (mW)','fontsize',14, 'interpreter', 'latex')
zlabel('P_{out}(\muW)')
view(-132.5,35)
colorbar
%%
figure
PinPoutTot(1,:)=Pin2(1)./PoutTot(1,:).*1000;
PinPoutTot(2,:)=Pin2(2)./PoutTot(2,:).*1000;
PinPoutTot(3,:)=Pin2(3)./PoutTot(3,:).*1000;
PinPoutTot(4,:)=Pin2(4)./PoutTot(4,:).*1000;
PinPoutTot(5,:)=Pin2(5)./PoutTot(5,:).*1000;
surf(X,Y,PinPoutTot)
shading interp
xlabel('$\phi^{(2)} (\rm ps^{2})$','fontsize',14,'interpreter', 'latex')
ylabel('Pin (mW)','fontsize',14, 'interpreter', 'latex')
zlabel('P_{in}/P_{out}')
view(-132.5,35)
colorbar
