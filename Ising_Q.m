clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 12; %size
dt = 1000;
T = 1000000;
t = 0:dt:T;
nt = length(t);

p = 2;
m = zeros(1,nt);

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

load('mem_N12_p2_No1.mat')
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
    H = H -Jij(1,j).*(H1 + H1');
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
        H = H -Jij(i,j).*(H1 + H1');
    end
end

% pos=1单独赋值
H2 = sigma_x;
for j = 2:N
    H2 = kron(H2,I2);
end
H = H - h.*H2;

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

temp = mem_con{1};
temp = round((temp+1)/2)+1;
spin0 = sigma_z(:,temp(1));
for i = 2:N
    spin0 = kron(spin0,sigma_z(:,temp(i)));
end

m(1) = 1;

% spin = spin0;
spin = gpuArray(spin0);
% H = gpuArray(H);
trans = expm(-1i*H*dt);
trans = gpuArray(trans);
for i = 2:nt
    spin = trans*spin;
    m(i) = spin'*spin0;
end
mm = abs(m);

figure;
plot(t,mm);



toc;