function [dtd,dtn,A_] = alertness_gen(sample_rate,WL,sample_mode,force_data)    
    
    %% Create the Alertness Model
    
    U=14.3; L=2.4;
    Sn = @(t,Sr,U) U-((U-Sr).*exp(-0.381*t));
    S  = @(t,Sa,L) (Sa-L).*exp(-0.0353*t)+L;
    C  = @(t) 2.52*cos((t-16.835)*pi/12);
    
    %% Create the time matrix
    N = size(WL,1);                              % numero de dias
    lim = reshape(WL',1,N*2);                    % carga linear
    tlim = cumsum(lim);                          % carga acumulada
    ptsPerHour = 100;
    %t = linspace(0,N*24,N*ptsPerHour*24); 
    
    %t = linspace(0,N*24,N*ptsPerHour*24); %resolution as 80 points per hour
    
    %% Determine the Alertness level
    % Circadian component
    %CP = C(t);
    
    % Homeostatic component
    Sr = zeros(N,1);
    Sa(1,:) = 14.3; %Sa(1,:) = N*2; initialize in a random value
    %Sa(1,:) = N*2;
    init = 0;
    for k=1:N

        td = round(linspace(0,WL(k,2),ptsPerHour*WL(k,2)),2); %+ init; %Day
        HCDW = S(td,Sa(k),L);
        tdr = td + init;%
        CPD = C(tdr);%
        
        Sr(k) = HCDW(end);
        tn = round(linspace(0,WL(k,1),ptsPerHour*WL(k,1)),2); %Night
        HCDS = Sn(tn,Sr(k),U);
        tnr = tn + tdr(end);%
        CPS = C(tnr);%
        
        Sa(k+1,:) = HCDS(end);    
        %aux(k,:) = [HCDW HCDS]; %A row for each day
        %aux1(k,:) = [HCDW zeros(size(HCDS))];   
        
        homeostatic{k} = [HCDW HCDS];%
        day_hom{k} = [HCDW zeros(size(HCDS))];%
        circadian{k} = [CPD CPS];%
        time_stamp{k} = [tdr tnr];%
        init = tnr(end);%
%         figure(12)
%         subplot(2,2,1); hold on;
%         plot(time_stamp{k},circadian{k},'k','LineWidth',1.2);
%         subplot(2,2,2); hold on;
%         plot(time_stamp{k},homeostatic{k},'k','LineWidth',1.2);
%         plot(time_stamp{k},day_hom{k},'r--','LineWidth',1.1);
%         subplot(2,2,3); hold on;
%         plot(time_stamp{k},circadian{k}+homeostatic{k},'k','LineWidth',1.2);
%         subplot(2,2,4); hold on;
%         scatter(k,init-(k-1)*24,'ro','LineWidth',1.2);
    end
    
%     subplot(2,2,4);hold on;
%     plot(1:k,ones(k,1)*24,'Color',[.8 .8 .8],'LineWidth',1.2);
    
    CP = cell2mat(circadian);%
    t = cell2mat(time_stamp);%
    
    %Initial values
    %SP = reshape(aux',1,N*ptsPerHour*24);
    %Day = reshape(aux1',1,N*ptsPerHour*24);

    SP = cell2mat(homeostatic);%
    Day = cell2mat(day_hom);%
    A = CP + SP;
    
    %% Debug figures
%     figure(13)
%     subplot(3,1,1)
%     plot(t,CP,'k','LineWidth',1.2);
%     subplot(3,1,2)
%     plot(t,SP,'k','LineWidth',1.2); hold on;
%     plot(t,Day,'r--','LineWidth',1.1);
%     subplot(3,1,3)
%     plot(t,A,'k','LineWidth',1.2);
    
    %% Determine the estimation data for nigth and day as windows
    
    p = ones(1,length(A));
    kw = find((SP-Day) == 0);
    p(kw) = 0;
    Aw = A(kw);     
    tw = t(kw);
    
    % Determine the day--night change moments
    chg = zeros(1,length(A));
    for i = 2 : length(A)
        if p(i) ~= p(i-1)
            chg(i) = 1;
        else
            chg(i) = 0;
        end
    end
    
    sampler = ptsPerHour;
    ptsPerHour = sample_rate;

    k_d = 1; k_n = 1;
    j_d = 1; j_n = 1;
    for i = 1 : length(A)
       if p(i) == 0 %IF day THEN
          A_.Ad{k_d}(j_d) = A(i);
          A_.td{k_d}(j_d) = t(i);
          j_d = j_d + 1;
       elseif p(i) == 1 %IF nigth THEN
          A_.An{k_n}(j_n) = A(i);
          A_.tn{k_n}(j_n) = t(i);
          j_n = j_n + 1;
       end

       if chg(i) == 1 && p(i) == 0
          k_n = k_n + 1; j_n = 1; %Reinitialize the night index
       elseif chg(i) == 1 && p(i) == 1
          k_d = k_d + 1; j_d = 1; % Reinitialize the day index
       end
    end
    
    for i = 1 : N
       
       %% Determine the points to reset at each part of the day
       
       % day indexes linear spaced
       ind_d = unique(round(linspace(1,length(A_.Ad{i}),...
                                              round(WL(i,2)*ptsPerHour))));    
       % night indexes linear spaced
       ind_n = unique(round(linspace(1,length(A_.An{i}),...
                                              round(WL(i,1)*ptsPerHour))));
       
       % if it is variable step approach
       if strcmp(sample_mode,'random')
          cruze = sampler*WL(i,2)/(length(ind_d)*ptsPerHour);
          cruze_rand = 2*rand(length(ind_d),1).*cruze - cruze;
          ind_d = sort(abs(unique(round(ind_d' + cruze_rand))));
          cruze_ran = 2*rand(length(ind_n),1).*cruze - cruze;
          ind_n = sort(abs(unique(round(ind_n' + cruze_ran))));
          ind_d = ind_d(ind_d <= length(A_.Ad{i})); 
          ind_n = ind_n(ind_n <= length(A_.An{i}));
          ind_d = ind_d(ind_d > 0);
          ind_n = ind_n(ind_n > 0);
       end                                     
       
       % if wants to force the window first point to the data
       if strcmp('first',force_data)
           % First data of window is included
           if ind_d(1) ~= 1
              ind_d = [1; ind_d];
           end
           if ind_n(1) ~= 1
              ind_n = [1; ind_n];
           end
       % if wants to force the window first and last points to the data
       elseif strcmp('first and last',force_data)
           % First data of window is included
           if ind_d(1) ~= 1
              ind_d = [1; ind_d];
           end
           if ind_n(1) ~= 1
              ind_n = [1; ind_n];
           end
           
           %Last window data is included
           if ind_d(end) ~= length(A_.Ad{i})
              ind_d = [ind_d; length(A_.Ad{i})];
           end
           if ind_n(end) ~= length(A_.An{i})
              ind_n = [ind_n; length(A_.An{i})];
           end
       % if wants to force the window last point to the data
       elseif strcmp('last',force_data)
           
           %Last window data is included
           if ind_d(end) ~= length(A_.Ad{i})
              ind_d = [ind_d; length(A_.Ad{i})];
           end
           if ind_n(end) ~= length(A_.An{i})
              ind_n = [ind_n; length(A_.An{i})];
           end
           
       end
       
       % Re-structure the data to easily estimate
       dtd.y{i} = A_.Ad{i}(ind_d)'; dtd.t{i} = A_.td{i}(ind_d)';
       dtd.initial(i) = A_.td{i}(1); dtd.final(i) = A_.td{i}(end);
       dtn.y{i} = A_.An{i}(ind_n)'; dtn.t{i} = A_.tn{i}(ind_n)';
       A_.Ad{i} = A_.Ad{i}'; A_.td{i} = A_.td{i}';
       A_.An{i} = A_.An{i}'; A_.tn{i} = A_.tn{i}';  
    end
    
    % Re-structure again the cells
    dtd.y = dtd.y'; dtd.t = dtd.t'; A_.Ad = A_.Ad'; A_.td = A_.td';
    dtn.y = dtn.y'; dtn.t = dtn.t'; A_.An = A_.An'; A_.tn = A_.tn';
    
    A_.A = A'; A_.t = t';
    
end