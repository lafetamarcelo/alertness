%% Sandro data-set
clear; close; clc;
load('sandro_all.mat');

%% 

date.time = datetime(kssApp.datetime);
date.diff = hours(diff(date.time));

ref_days = floor(days(date.time- min(date.time)));

f_day = day(min(date.time));
f_month = month(min(date.time));
date_sec = string(f_day) + '/' + string(f_month) + '/2018 04:00:00';
init = datetime(date_sec, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');

d_days = unique(floor(days((date.time-min(date.time)))));
ind = min(d_days):max(d_days);
graph.time = init + hours(24.*ind)';

%Figure plot KSS as sleepiness
figure(1)
hax=axes; 
scatter(date.time, kssApp.value, 'r', 'LineWidth', 1.4); 
hold on;
for i = 1 : length(graph.time)
    line([graph.time graph.time], get(hax,'YLim'),...
            'Color', [.6 .6 .6], 'LineWidth', 1.2);
end
%plot(date.time, flip(kss), 'r--', 'LineWidth', 1.2);
hold off;

alertness = (10 - kssApp.value).*(13/9) + 1;

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

aux = kssApp.value; k = 1; w = 1;
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

%% simple pre treatment

% Filter out specific days
filt_days = 15; k = 1;
for i = 1 : length(KSS.time)
   if not(ismember(i,filt_days))
       KSS.time_filt{k} = KSS.time{i};
       KSS.kss_filt{k} = KSS.kss{i};
       KSS.alert_filt{k} = KSS.alert{i};
       k = k + 1;
   end
end

% Transform the time data into hours...
time_ref = KSS.time_filt{1}(2);
for i = 1 : length(KSS.time_filt)
    KSS.y{i,1} = KSS.alert_filt{i}(2:end)';
    KSS.t{i,1} = hours(KSS.time_filt{i}(2:end)' - time_ref);
end

% Ploting the resulted manipulations
figure(2); hold on;
for i = 1 : length(KSS.time_filt)
    plot(KSS.time_filt{i}, KSS.alert_filt{i}, 'k.-', 'LineWidth', 1.2);
    plot(hours(KSS.t{i}) + time_ref, KSS.y{i} + 1, 'k--', 'LineWidth', 1.2)
end
hold off;

%% save only interesting data

dte.t = KSS.t;
dte.y = KSS.y;
for i = 1 : length(dte.t)
   dte.init(i) = dte.t{i}(1); 
   dte.final(i) = dte.t{i}(end);
end
save('regres_sandro.mat', 'dte')