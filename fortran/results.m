clear
clc

P = load('parameters.txt');
M = load('magnetization.txt');
E = load('energy.txt');
S = load('state_evolution.txt');
L = P(1);
N = P(2);
f = P(3);

Si(1:L,1:L,1:N/f) = 0;
n = 0;
for it = 1:f:N
    for i = 1:L
        for j = 1:L
            n = n+1;
            Si(i,j,it) = S(n); 
        end
    end
end

figure(1)
for i = 1:f:N
    imagesc(Si(:,:,i))
    title(sprintf('Spin Configuration at step %d',i))
    xlabel('x')
    ylabel('y')
    pause(0.01)
end

figure(2)
plot(1:f:N, M, 'LineWidth', 1.5)
title('Magnetization vs Monte Carlo Steps');
xlabel('Monte Carlo Steps');
ylabel('Magnetization');

figure(3)
plot(1:f:N, E, 'LineWidth', 1.5);
title('Energy vs Monte Carlo Steps');
xlabel('Monte Carlo Steps');
ylabel('Energy');