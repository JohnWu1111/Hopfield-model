clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 8; %size
dt = 1;
T = 100;
t = 0:dt:T;
nt = length(t);

p = 3;

h = 0.1;

% Cv = zeros(1,round((Tmax-Tmin)/Tstep+1));

% parameters of memeory

mem_con = cell(p,1);
Jij = zeros(N,N);
% for k = 1:p
%     mem_con{k} = round(rand(1,N))*2-1;
%     temp = mem_con{k};
%     for i = 1:N
%         for j = i+1:N
%             Jij(i,j) = Jij(i,j) + temp(i)*temp(j);
%         end
%     end
% end
mem_con{1} = ones(1,N);
mem_con{2} = [ones(1,6) -ones(1,2)];
mem_con{3} = [ones(1,3) -ones(1,5)];
for k = 1:p
    temp = mem_con{k};
    for i = 1:N
        for j = i+1:N
            Jij(i,j) = Jij(i,j) + temp(i)*temp(j);
        end
    end
end
Jij = Jij/N;

% load('mem_N12_p3_No1.mat')
% overlap12 = mean(mem_con{1}.*mem_con{2})
% overlap13 = mean(mem_con{1}.*mem_con{3})
% overlap23 = mean(mem_con{2}.*mem_con{3})

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
    H = H -Jij(1,j).*H1/4;
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
        H = H -Jij(i,j).*H1/4;
    end
end

% pos=1单独赋值
H2 = sigma_x;
for j = 2:N
    H2 = kron(H2,I2);
end
H = H + h.*H2/2;

for i = 2:N
    H2 = I2;
    for j = 2:i-1
        H2 = kron(H2,I2);
    end
    H2 = kron(H2,sigma_x);
    for j = i+1:N
        H2 = kron(H2,I2);
    end
    H = H + h.*H2/2;
end

temp1 = mem_con{1};
temp1 = round((temp1+1)/2)+1;
spin1 = sigma_x(:,temp1(1));
temp2 = mem_con{2};
temp2 = round((temp2+1)/2)+1;
spin2 = sigma_x(:,temp2(1));
temp3 = mem_con{3};
temp3 = round((temp3+1)/2)+1;
spin3 = sigma_x(:,temp3(1));
for i = 2:N
    spin1 = kron(spin1,sigma_x(:,temp1(i)));
    spin2 = kron(spin2,sigma_x(:,temp2(i)));
    spin3 = kron(spin3,sigma_x(:,temp3(i)));
end

% m1 = zeros(1,nt);
% m1(1) = 1;
% m2 = zeros(1,nt);
% m2(1) = spin2'*spin1;
% m3 = zeros(1,nt);
% m3(1) = spin3'*spin1;

% spin = gpuArray(spin1);
% trans = expm(-1i*H*dt);
% trans = gpuArray(trans);
% for i = 2:nt
%     spin = trans*spin;
%     m1(i) = spin'*spin1;
%     m2(i) = spin'*spin2;
%     m3(i) = spin'*spin3;
% end

% m1 = abs(m1);
% m2 = abs(m2);
% m3 = abs(m3);

H = full(H);
H = gpuArray(H);
[V,D] = eig(H);
e = diag(D);
spin0 = V'*spin1;
trans = exp(-1i*e*t);
spin = spin0.*trans;
spint = V*spin;
m1 = abs(spin1'*spint);
m2 = abs(spin2'*spint);
m3 = abs(spin3'*spint);



figure;
plot(t,m1,t,m2,t,m3);



toc;