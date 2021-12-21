clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 12; %size
dt = 1;
T = 1000;
t = 0:dt:T;
nt = length(t);
J = 1;

h = 0.3;

% sigma_x = sparse([0,1;1,0])./2;
% sigma_z = sparse([1,0;0,-1])./2;
sigma_x = sparse([0,1;1,0]);
sigma_z = sparse([1,0;0,-1]);
I2 = sparse(eye(2));

H = sparse(2^N,2^N);

% pos=1单独赋值
for j = 2:N
    H1 = sigma_z;
    for k = 2:j-1
        H1 = kron(H1,I2);
    end
    H1 = kron(H1,sigma_z);
    for k = j+1:N
        H1 = kron(H1,I2);
    end
    H = H -H1*J/N;
end

for i = 2:N-1    
    for j = i+1:N
        H1 = I2;
        for k = 2:i-1
            H1 = kron(H1,I2);
        end
        H1 = kron(H1,sigma_z);
        for k = i+1:j-1
            H1 = kron(H1,I2);
        end
        H1 = kron(H1,sigma_z);
        for k = j+1:N
            H1 = kron(H1,I2);
        end
        H = H -H1*J/N;
    end
end

% pos=1单独赋值
H2 = sigma_x;
for j = 2:N
    H2 = kron(H2,I2);
end
H = H + h.*H2;

for i = 2:N
    H2 = I2;
    for j = 2:i-1
        H2 = kron(H2,I2);
    end
    H2 = kron(H2,sigma_x);
    for j = i+1:N
        H2 = kron(H2,I2);
    end
    H = H + h.*H2;
end

spin1 = zeros(2^N,1);
spin1a = zeros(2^N,1);
spin1(1) = 1;
spin1a(end) = 1;

m1 = zeros(1,nt);
m1a = zeros(1,nt);
m1(1) = 1;
m1a(1) = 0;

% spin = spin0;
spin = gpuArray(spin1);
% H = gpuArray(H);
trans = expm(-1i*H*dt);
trans = gpuArray(trans);
for i = 2:nt
    spin = trans*spin;
    m1(i) = spin'*spin1;
    m1a(i) = spin'*spin1a;
end
m1 = abs(m1);
m1a = abs(m1a);

figure;
plot(t,m1,t,m1a);
legend('m1','m1a');

toc;