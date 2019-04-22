%% Read the data-set
clear; close; clc;
importKSS;

%% Data preparation for KSS and PVT

time_sampled = string(data) + ' ' + string(hora);
date.time = datetime(flip(time_sampled), ...
                        'InputFormat','dd/MM/yyyy HH:mm:ss');
date.diff = hours(diff(date.time));

f_day = day(min(date.time));
f_month = month(min(date.time));
date_sec = string(f_day) + '/' + string(f_month) + '/2019 02:00:00';
init = datetime(date_sec, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');

d_days = unique(floor(days((date.time-min(date.time)))));
ind = min(d_days):max(d_days);
graph.time = init + hours(24.*ind)';

%Figure plot KSS as sleepiness
figure(1)
hax=axes; 
scatter(date.time, flip(kss), 'r', 'LineWidth', 1.4); 
hold on;
for i = 1 : length(graph.time)
    line([graph.time graph.time], get(hax,'YLim'),...
            'Color', [.6 .6 .6], 'LineWidth', 1.2);
end
%plot(date.time, flip(kss), 'r--', 'LineWidth', 1.2);
hold off;

%Figure plot KSS as alertness

alertness = (10 - flip(kss)).*(13/9) + 1;

figure(2)
hax = axes; 
scatter(date.time, alertness, 'r', 'LineWidth', 1.4); 
hold on;
for i = 1 : length(graph.time)
    line([graph.time graph.time], get(hax,'YLim'),...
            'Color', [.6 .6 .6], 'LineWidth', 1.2);
end
%plot(date.time, alertness, 'r--', 'LineWidth', 1.2);
hold off;

% Determine the KSS alertness data pack

aux = flip(kss); k = 1; w = 1;
flag = 0;
for i = 1 : length(graph.time)
    for j = 1 : length(date.time)
        if i ~= length(graph.time)
            if (date.time(j) >= graph.time(i)) && (date.time(j) <= graph.time(i+1))
                KSS.kss{w}(k) = aux(j);
                KSS.alert{w}(k) = (10 - aux(j))*(13/9) + 1;
                KSS.time{w}(k) = date.time(j);
                k = k + 1; flag = 1;
            end
        elseif date.time(j) >= graph.time(i)
            KSS.kss{w}(k) = aux(j);
            KSS.alert{w}(k) = (10 - aux(j))*(13/9) + 1;
            KSS.time{w}(k) = date.time(j);
            k = k + 1; flag = 1;
        end
    end
    k = 1;
    
    if flag == 1
       w = w + 1;
       flag = 0;
    end
end

% Certify the data manipulation
for i = 1 : length(KSS.time)
    figure(1); hold on;
    scatter(KSS.time{i}, KSS.kss{i}, 'x', 'LineWidth', 1.2); hold off;
    figure(2); hold on;
    scatter(KSS.time{i}, KSS.alert{i}, 'x', 'LineWidth', 1.2); hold off;
end

% Plotting the difference between windows
figure(3)
hax = axes; 
plot(date.time(2:end), date.diff, 'k-', 'LineWidth', 1.4); hold on;
for i = 1 : length(graph.time)
    line([graph.time graph.time], get(hax,'YLim'),...
                     'Color', [.6 .6 .6], 'LineWidth', 1.2);
end
hold off;

%% Minor data treatments

% Filter out specific days
filt_days = [1,2]; k = 1;
for i = 1 : length(KSS.time)
   if not(ismember(i,filt_days))
       KSS.time_filt{k} = KSS.time{i};
       KSS.kss_filt{k} = KSS.kss{i};
       KSS.alert_filt{k} = KSS.alert{i};
       k = k + 1;
   end
end

% Transform the time data into hours...
time_ref = KSS.time_filt{1}(1);
for i = 1 : length(KSS.time_filt)
    KSS.y{i,1} = KSS.alert_filt{i}';
    KSS.t{i,1} = hours(KSS.time_filt{i}' - time_ref);
end

% Ploting the resulted manipulations
figure(2); hold on;
for i = 1 : length(KSS.time_filt)
    plot(KSS.time_filt{i}, KSS.alert_filt{i}, 'k.-', 'LineWidth', 1.2);
    plot(hours(KSS.t{i}) + time_ref, KSS.y{i} + 1, 'k--', 'LineWidth', 1.2)
end
hold off;

%% Determine possible resampling of the data

desiredFs = 1;

% Low pass filtering design
p = 1;
q = 4;

% ensure an odd length filter
n = 10*q+1;

% use .25 of Nyquist range of desired sample rate
cutoffRatio = .35;

% construct lowpass filter 
lpFilt = p * fir1(n, cutoffRatio * 1/q);


figure(4); hold on;
for i = 1 : length(KSS.t)
    x = KSS.y{i}; t = KSS.t{i} - KSS.t{i}(1);
    a(1) = (x(end)-x(1)) / (t(end)-t(1));
    a(2) = x(1);
    
    xdetrend = x - polyval(a,t);
   
    [ydetrend,ty] = resample(xdetrend,t,desiredFs,p,q,lpFilt);
    y = ydetrend + polyval(a,ty); 
    
    KSS.re_y{i,1} = y;
    KSS.re_t{i,1} = ty + KSS.t{i}(1);
    scatter(KSS.t{i}, KSS.y{i}, 'ko', 'LineWidth', 1.4);
    plot(KSS.re_t{i}, KSS.re_y{i}, 'r-', 'LineWidth', 1.2);
    plot(KSS.re_t{i}, KSS.re_y{i}, 'r.', 'LineWidth', 1.2);
end

%% save only interesting data

dte.t = KSS.t;
dtr.realTime = KSS.time_filt;
dtr.timeRef = timeRef;
dte.y = KSS.y;
for i = 1 : length(dte.t)
   dte.init(i) = dte.t{i}(1); 
   dte.final(i) = dte.t{i}(end);
end
save('regres_02.mat', 'dte');

%% plot the reference for determining the initial and final for each day

figure(5);
subplot(2,1,1);
scatter(cell2mat(KSS.t), cell2mat(KSS.y), 'ko', 'LineWidth', 1.4);
subplot(2,1,2); hold on;
for i = 1 : length(KSS.time_filt)
    scatter(KSS.time_filt{i}, KSS.alert_filt{i}, 'ko', 'LineWidth', 1.4);
end
gInd = find(graph.time >= KSS.time_filt{1}(1));
line([graph.time(gInd) graph.time(gInd)], [14 4],...
                     'Color', [.6 .6 .6], 'LineWidth', 1.2);
%% Determine with particular different initials and finals

dte.t = KSS.t; dte.y = KSS.y;
dtr.realTime = KSS.time_filt;
dtr.timeRef = time_ref;
removeInit = [0 4 4 5];
addFinal = [0 0 0 0];

for i =  1 : length(dte.t)
     dte.init(i) = dte.t{i}(1) - removeInit(i);
     dte.final(i) = dte.t{i}(end) + addFinal(i);
     figure(5); hold on;
     subplot(2,1,1); hold on;
     line([dte.init(i) dte.init(i)], [14 4], 'Color', [1 0 0], 'LineWidth', 1.4);
     line([dte.final(i) dte.final(i)], [14 4], 'Color', [0 0 1], 'LineWidth', 1.4);
end
%%
save('regres_02.mat', 'dte');