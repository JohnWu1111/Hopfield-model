clear;
% close all;
clc;
format long
tic;

% Definition of parameters
n = 20; %size
N = n^2;
tmax = 1e7; % eventnumber

T = 0;

skip = 0;

% Cv = zeros(1,round((Tmax-Tmin)/Tstep+1));

% parameters of memeory
p = 2;
mem_con = cell(p,1);
coeff_save = zeros(N,N);
coeff_save2 = cell(n,n);
for i = 1:p
    mem_con{i} = round(rand(n))*2-1;
    coeff_save = coeff_save + kron(mem_con{i},mem_con{i});
end

for i = 1:n
    for j = 1:n
        coeff_save2{i,j} = coeff_save((i-1)*n+1:i*n,(j-1)*n+1:j*n);
    end
end

spin = mem_con{1};

Ett = 0;
% Activation
for s = 1:skip
    r = ceil(rand(1,2)*n);
    coeff_it = coeff_save2{r(1),r(2)};
    dEt = (mean(coeff_it.*spin,'all')-coeff_it(r(1),r(2))*spin(r(1),r(2))/N)*spin(r(1),r(2));
    
    if rand < exp(-dEt/T)
        spin(r(1),r(2)) = -spin(r(1),r(2));
        Ett = Ett + dEt;
    end
end

spintotal = zeros(p,tmax);
%     Etotal = zeros(1,tmax);

for i = 1:tmax %模拟次数循环
    r = ceil(rand(1,2)*n);
    coeff_it = coeff_save2{r(1),r(2)};
    dEt = (mean(coeff_it,'all')-coeff_it(r(1),r(2))/N)*spin(r(1),r(2));
    
    if rand < exp(-dEt/T) %判断是否作出状态改变
        spin(r(1),r(2)) = -spin(r(1),r(2));
        Ett = Ett + dEt;
    end
    
    %         Etotal(i) = Ett;
    for j = 1:p
        temp = mem_con{j};
        spintotal(j,i) = mean(spin.*temp,'all')^2;
    end
end

m2 = sqrt(mean(spintotal,2))
map = zeros(n+1,n+1);

for k = 1:p
    figure;
    map(1:n,1:n) = mem_con{k};
    pcolor(double(map));
    colormap(gray(2));
end
figure;
map(1:n,1:n) = spin;
pcolor(double(map));
colormap(gray(2));

%     Cv(j) = var(Etotal)/T(j)^2;


toc;