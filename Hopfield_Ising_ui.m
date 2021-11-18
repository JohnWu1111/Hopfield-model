clear;
% close all;
clc;
format long
tic;

% Definition of parameters
n = 50; %size

tmax = 1e7; % eventnumber
tstep = 1e4;
nt = round(tmax/tstep);
t = 1:nt;

Tmax = 3; %temperature
Tmin = 0.5;
Tstep = 0.05;
T = Tmin:Tstep:Tmax;

skip = 1e7;

spin = round(rand(n))*2-1;
spint = zeros(n+1);
spint(1:n,1:n) = spin;
Et = zeros(1,nt);

Cv = zeros(1,round((Tmax-Tmin)/Tstep+1));

% parameters of memeory
N = 2;
mem_con = cell(2,1);
for i = 1:N
    mem_con{i} = round(rand(n))*2-1;
end

% Initialization of figures
figure;
subplot(2,2,1)
gspin = pcolor(double(spint));
colormap(gray(2));

subplot(2,2,2)
gEt = plot(t,Et);
title('E')
xlabel('eventnumber')
ylabel('E')

subplot(2,2,3)
gCv = plot(T,Cv,'-o');
title('Cv')
xlabel('T')
ylabel('Cv')

% Activation
for s=1:skip
    r = ceil(rand(1,2)*n);
    su = spin(mod(r(1)-2,n)+1,r(2));
    sd = spin(mod(r(1),n)+1,r(2));
    sl = spin(r(1),mod(r(2)-2,n)+1);
    sr = spin(r(1),mod(r(2),n)+1);
    it = spin(r(1),r(2));
    dEt = 0;
    for k = 1:N
        temp = mem_con{k};
        mu = temp(mod(r(1)-2,n)+1,r(2));
        md = temp(mod(r(1),n)+1,r(2));
        ml = temp(r(1),mod(r(2)-2,n)+1);
        mr = temp(r(1),mod(r(2),n)+1);
        mm = temp(r(1),r(2));
        dEt = dEt + 2*it*su*mm*mu;
        dEt = dEt + 2*it*sd*mm*md;
        dEt = dEt + 2*it*sl*mm*ml;
        dEt = dEt + 2*it*sr*mm*mr;
    end
    dEt = dEt/N;
    %     dEt = spin(r(1),r(2))*(su+sd+sl+sr)*2;
    if rand < exp(-dEt/Tmin)
        spin(r(1),r(2)) = -spin(r(1),r(2));
    end
end

% Temperature loop
Tii = 1;
for T = Tmin:Tstep:Tmax
    spintotal = zeros(1,tmax);
    Etotal = zeros(1,tmax);
    
    %     idx = repmat({':'},ndims(spin),1);
    %     idx{1} = [n 1:n-1];
    %     nspin = spin(idx{:});
    %     idx = repmat({':'},ndims(spin),1);
    %     idx{1} = [2:n 1];
    %     nspin = nspin+spin(idx{:});
    %
    %     idx = repmat({':'},ndims(spin),1);
    %     idx{2} = [n 1:n-1];
    %     nspin = nspin+spin(idx{:});
    %     idx = repmat({':'},ndims(spin),1);
    %     idx{2} = [2:n 1];
    %     nspin = nspin+spin(idx{:});
    %     Ett = sum(sum(spin.*nspin));
    Ett = 0;
    
    j = 1;
    for  i = 1:tmax %模拟次数循环
        r = ceil(rand(1,2)*n);
        su = spin(mod(r(1)-2,n)+1,r(2));
        sd = spin(mod(r(1),n)+1,r(2));
        sl = spin(r(1),mod(r(2)-2,n)+1);
        sr = spin(r(1),mod(r(2),n)+1);
        %         dEt = spin(r(1),r(2))*(su+sd+sl+sr)*2;
        it = spin(r(1),r(2));
        dEt = 0;
        for k = 1:N
            temp = mem_con{k};
            mu = temp(mod(r(1)-2,n)+1,r(2));
            md = temp(mod(r(1),n)+1,r(2));
            ml = temp(r(1),mod(r(2)-2,n)+1);
            mr = temp(r(1),mod(r(2),n)+1);
            mm = temp(r(1),r(2));
            dEt = dEt + 2*it*su*mm*mu;
            dEt = dEt + 2*it*sd*mm*md;
            dEt = dEt + 2*it*sl*mm*ml;
            dEt = dEt + 2*it*sr*mm*mr;
        end
        dEt = dEt/N;
        
        if rand < exp(-dEt/T) %判断是否作出状态改变
            spin(r(1),r(2)) = -spin(r(1),r(2));
            Ett = Ett + dEt;
        end
        
        Etotal(i) = Ett;
        %         spintotal(i) = mean(mean(spin));
        
        if mod(i,tstep) == 0
            Et(j) = Ett;
            spint(1:n,1:n) = spin;
            set(gspin,'Cdata',double(spint));
            set(gEt,'Ydata',double(Et));
            drawnow
            j = j+1;
        end
    end
    Cv(Tii) = var(Etotal)/T^2;
    Tii = Tii+1;
    set(gCv,'Ydata',double(Cv));
    drawnow
    
    %     if (T<2) && (Mtime(Tii-1)<0.8) %排除不应该出现的解
    %         error('模拟不理想，请再运行一次')
    %     end
    
end

toc;