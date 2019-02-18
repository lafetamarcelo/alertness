function windows = autoGen(days)
   
    wake = round(sqrt(2/3)*randn(days,1) + 16,1);
    sleep = 24 - wake;
    windows = [sleep, wake];
    
    %% Debug figure
%     figure(1)
%     subplot(3,1,1)
%     title({['E(x) =',num2str(mean(wake))]})
%     plot(1:days,ones(days,1)*16,'Color',[.8 .8 .8],...
%         'LineWidth',1.2);hold on;
%     plot(1:days,wake,'r','LineWidth',1.5);
%     subplot(3,1,2)
%     title(['E(x) =',num2str(mean(sleep))])
%     plot(1:days,ones(days,1)*8,'Color',[.8 .8 .8],...
%         'LineWidth',1.2);hold on;
%     plot(1:days,sleep,'r','LineWidth',1.5);
%     subplot(3,1,3)
%     plot(1:days,ones(days,1)*24,'Color',[.8 .8 .8],...
%         'LineWidth',1.2);hold on;
%     plot(1:days,wake+sleep,'r','LineWidth',1.5);
%     
%     disp(['E(x) =',num2str(mean(wake))])
%     disp(['E(x) =',num2str(mean(sleep))])
%     disp(['E(x) =',num2str(mean(wake+sleep))])
    
end