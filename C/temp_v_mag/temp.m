clear
clc

data = load("data.txt");
N = 200*200;

figure(1)
plot(abs(data(:,2)),abs(data(:,1)/N),LineWidth=1.5)