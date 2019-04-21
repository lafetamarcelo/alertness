function windows = autoGen(days)
   
    wake = round(sqrt(2/3)*randn(days,1) + 16,1);
    sleep = 24 - wake;
    windows = [sleep, wake];    
end