clear;
% close all;
clc;
format long
tic;

% Definition of parameters
n = 5; %size
N = n^2;
tmax = 1e3; % eventnumber
tstep = 1;
nt = round(tmax/tstep);
t=1:nt;

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

% load('mem_N16_p10_No1.mat');

spin = mem_con{1};

m2t = zeros(p,nt);
Et = zeros(1,nt);
figure;

subplot(2,2,1)
gM1 = plot(t,m2t(1,:));

subplot(2,2,2)
gM2 = plot(t,m2t(2,:));

% subplot(2,2,3)
% gM3 = plot(t,m2t(3,:));

subplot(2,2,4)
gEt = plot(t,Et);

Ett = 0;
ii = 1;

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
Etotal = zeros(1,tmax);

for i = 1:tmax %模拟次数循环
    r = ceil(rand(1,2)*n);
    coeff_it = coeff_save2{r(1),r(2)};
    dEt = (mean(coeff_it.*spin,'all')-coeff_it(r(1),r(2))*spin(r(1),r(2))/N)*spin(r(1),r(2));
    
    if rand < exp(-dEt/T) %判断是否作出状态改变
        spin(r(1),r(2)) = -spin(r(1),r(2));
        Ett = Ett + dEt;
    end
    
    Etotal(i) = Ett;
    for j = 1:p
        temp = mem_con{j};
        spintotal(j,i) = mean(spin.*temp,'all')^2;
    end
    
    if mod(i-1,tstep)==0
        m2t(:,ii) = sqrt(spintotal(:,i));
        Et(ii) = Ett;
        ii = ii +1;
        set(gM1,'Ydata',double(m2t(1,:)));
        set(gM2,'Ydata',double(m2t(2,:)));
%         set(gM3,'Ydata',double(m2t(3,:)));
        set(gEt,'Ydata',double(Et));
        drawnow
    end
end

m2 = sqrt(mean(spintotal,2))

% for k = 1:N
%     figure;
%     temp = mem_con{k};
%     pcolor(double(temp));
%     colormap(gray(2));
% end
% figure;
% pcolor(double(spin));
% colormap(gray(2));

% Cv = var(Etotal)/T


toc;