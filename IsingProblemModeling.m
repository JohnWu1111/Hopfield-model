clear all; %格式化及格式配置
close all;
clc;
format long

n=50; %初始化变量
tmax=1e6; 
Tmax=3; 
Tmin=1; 
tstep=1e4; 
Tstep=1/50; 
skip=1e6; 
nt=round(tmax/tstep); 
t=1:nt; 
T=Tmin:Tstep:Tmax;
Tii=1;
spin = round(rand(n))*2-1;
spint = zeros(n+1);
spint(1:n,1:n)=spin;
Mt=zeros(1,nt);
Et=zeros(1,nt);
Mtime=zeros(1,round((Tmax-Tmin)/Tstep+1));
Cv=zeros(1,round((Tmax-Tmin)/Tstep+1));
Kai=zeros(1,round((Tmax-Tmin)/Tstep+1));

figure; %图表初始化

subplot(3,2,1)
gspin=pcolor(double(spint));
colormap(gray(2));

subplot(3,2,2)
gMt=plot(t,Mt);
title('M')
xlabel('eventnumber')
ylabel('M')

subplot(3,2,3)
gEt=plot(t,Et);
title('E')
xlabel('eventnumber')
ylabel('E')

subplot(3,2,4)
gMTime=plot(T,Mtime,'-o');
title('<M>')
xlabel('T')
ylabel('<M>')

subplot(3,2,5)
gCv=plot(T,Cv,'-o');
title('Cv')
xlabel('T')
ylabel('Cv')

subplot(3,2,6)
gKai=plot(T,Kai,'-o');
title('X')
xlabel('T')
ylabel('X')

 for s=1:skip  %初始平衡过程
    r=ceil(rand(1,2)*n);
    su=spin(mod(r(1)-2,n)+1,r(2));
    sd=spin(mod(r(1),n)+1,r(2));
    sl=spin(r(1),mod(r(2)-2,n)+1);
    sr=spin(r(1),mod(r(2),n)+1);
    dEt=spin(r(1),r(2))*(su+sd+sl+sr)*2;
    if rand<exp(-dEt/Tmin)
        spin(r(1),r(2))=-spin(r(1),r(2));
    end
 end
 
for T=Tmin:Tstep:Tmax %温度循环
spintotal=zeros(1,tmax);
Etotal=zeros(1,tmax);

idx=repmat({':'},ndims(spin),1);
idx{1}=[n 1:n-1];
nspin=spin(idx{:});
idx=repmat({':'},ndims(spin),1);
idx{1}=[2:n 1];
nspin=nspin+spin(idx{:});

idx=repmat({':'},ndims(spin),1);
idx{2}=[n 1:n-1];
nspin=nspin+spin(idx{:});
idx=repmat({':'},ndims(spin),1);
idx{2}=[2:n 1];
nspin=nspin+spin(idx{:});
Ett=sum(sum(spin.*nspin));

j=1;
for  i=1:tmax %模拟次数循环
    r=ceil(rand(1,2)*n);
    su=spin(mod(r(1)-2,n)+1,r(2));
    sd=spin(mod(r(1),n)+1,r(2));
    sl=spin(r(1),mod(r(2)-2,n)+1);
    sr=spin(r(1),mod(r(2),n)+1);
    dEt=spin(r(1),r(2))*(su+sd+sl+sr)*2;
    
    if rand<exp(-dEt/T) %判断是否作出状态改变
        spin(r(1),r(2))=-spin(r(1),r(2));
        Ett=Ett+dEt;
    end
    
    Etotal(i)=Ett;
    spintotal(i)=mean(mean(spin));
    
    if mod(i,tstep)==0 
        Mt(j)=sum(sum(spin));
        Et(j)=Ett;
        spint(1:n,1:n)=spin;
        set(gspin,'Cdata',double(spint));
        set(gMt,'Ydata',double(Mt));
        set(gEt,'Ydata',double(Et));
        drawnow
        j=j+1;
    end
end
Mtime(Tii)=abs(mean(spintotal));
Cv(Tii)=var(Etotal)/T^2;
Kai(Tii)=var(spintotal)/T;
Tii=Tii+1;
set(gMTime,'Ydata',double(Mtime));
set(gCv,'Ydata',double(Cv));
set(gKai,'Ydata',double(Kai));
drawnow

if (T<2) && (Mtime(Tii-1)<0.8) %排除不应该出现的解
   error('模拟不理想，请再运行一次') 
end

end
