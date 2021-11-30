clear;
% close all;
clc;
format long
tic;

% global dt H

maxNumCompThreads(1);

% Definition of parameters
N = 16; %size
dt = 0.1;
T = 10;
t = 0:dt:T;
nt = length(t);

p = 2;

h = 0.1;

% Cv = zeros(1,round((Tmax-Tmin)/Tstep+1));

% parameters of memeory

% mem_con = cell(p,1);
% Jij = zeros(N,N);
% for k = 1:p
%     mem_con{k} = round(rand(1,N))*2-1;
%     temp = mem_con{k};
%     for i = 1:N
%         for j = i+1:N
%             Jij(i,j) = Jij(i,j) + temp(i)*temp(j);
%         end
%     end
% end
% Jij = Jij/N;

load('mem_N16_p2_No1.mat')
% overlap = mean(mem_con{1}.*mem_con{2})

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
    H = H -Jij(1,j).*H1;
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
        H = H -Jij(i,j).*H1;
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

temp1 = mem_con{1};
temp1 = round((temp1+1)/2)+1;
spin1 = sigma_z(:,temp1(1));
temp2 = mem_con{2};
temp2 = round((temp2+1)/2)+1;
spin2 = sigma_z(:,temp2(1));
for i = 2:N
    spin1 = kron(spin1,sigma_z(:,temp1(i)));
    spin2 = kron(spin2,sigma_z(:,temp2(i)));
end

m1 = zeros(1,nt);
m1(1) = 1;
m2 = zeros(1,nt);
m2(1) = spin2'*spin1;

spin = full(spin1);
% spin = gpuArray(spin0);
% H = gpuArray(H);
% trans = expm(-1i*H*dt);
% trans = gpuArray(trans);
for i = 2:nt
%     spin = runge(spin);
    c1 = -1i*H*spin;
    c2 = -1i*H*(spin+c1.*dt/2);
    c3 = -1i*H*(spin+c2.*dt/2);
    c4 = -1i*H*(spin+c3.*dt);
    spin = spin + dt*(c1+2*c2+2*c3+c4)/6;
    m1(i) = spin'*spin1;
    m2(i) = spin'*spin2;
end

m1 = abs(m1);
m2 = abs(m2);

figure;
plot(t,m1,t,m2);

toc;

% function y = f(x)
%     global H
%     y = -1i*H*x;
% end
% 
% function y = runge(y0) 
%     global dt
%     y = y0;
%     
%     c1 = f(y);
%     c2 = f(y+c1.*dt/2);
%     c3 = f(y+c2.*dt/2);
%     c4 = f(y+c3.*dt);
%     y = y + dt*(c1+2*c2+2*c3+c4)/6;
% end
% 
% function y = runge(y0) 
%     global dt H
%     y = y0;
%     
%     c1 = -1i*H*y;
%     c2 = -1i*H*(y+c1.*dt/2);
%     c3 = -1i*H*(y+c2.*dt/2);
%     c4 = -1i*H*(y+c3.*dt);
%     y = y + dt*(c1+2*c2+2*c3+c4)/6;
% end