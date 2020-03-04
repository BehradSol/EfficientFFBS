clc
clear 
close all



T = 50;
p = 1;
Nx = 5;
Ny = 3;
Ne = 1;

A = cell(p,1);
for i = 1 : p
    A{i} = zeros(Nx,Nx);
end
A{1}(3,1) = -0.2;
A{1}(1,3) = 0.4;

A{1}(1,1) =  -0.75;
% A{2}(3,3) = -0.5;

B = zeros(Nx, Ne);
B(1,1) = 1;
e = rand(Ne , T);

qd = (1e-5)*ones(1,Nx);
qd(1) = 6;
qd(3) = 4;
Q = diag(qd);

R = 0.5*eye(Ny);

C = 2*rand(Ny,Nx)-1;

x = zeros(Nx,T+p);

for i = 1 : T
    x(:,i+p) = x(:,i+p) + B*e(:,i) + mvnrnd(zeros(1,Nx),Q)';
    for j = 1 : p
        x(:,i+p) = x(:,i+p) + A{j}*x(:,i+p-j);
    end
end

y = C*x(:,p+1 : T+p) + mvnrnd(zeros(1,Ny),R,T)';




[m,  Cov]   = EFBS(y, e, A, Q, B, C, R);
[m2, Cov2] = Filtering(y, e, A, Q, B, C, R);

plot(1:T , x(1,p+1:T+p) , 'black' , 'LineWidth' , 1.5)
hold on
plot(1:T , m2(1,:) , 'blue' , 'LineWidth' , 2.5)
plot(1:T , m(1,:) , 'red', 'LineWidth' , 1.5)

ylabel('Signal amplitude')
xlabel('Time index')

legend('Gound truth','Estimation via conventional filtering','Estimation via EFBS')
grid on
xlim([1 T])


figure

plot(1:T , x(1,p+1:T+p) , 'black' , 'LineWidth' , 1.5)
hold on
plot(1:T , m(1,:) , 'red' , 'LineWidth' , 1.5)
plot(1:T , m(1,:) + 3*Cov(1,1) , '--blue', 'LineWidth' , 1)
plot(1:T , m(1,:) - 3*Cov(1,1) , '--blue', 'LineWidth' , 1)


ylabel('Signal amplitude')
xlabel('Time index')

legend('Gound truth','Estimation via EFBS','95% Quantiles')
grid on
xlim([1 T])