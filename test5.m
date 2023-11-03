clear;

% figure;
% load('C:\Users\Administrator\Desktop\result_total_nondeg_p4_h0.1_res100_num1000.mat')
% count = 1;
% for i = 1:num
%     len = length(q_store{4,i});
%     a = round(0.45*len);
%     b = round(0.55*len);
%     temp1(count:count+b-a) = q_store{4,i}(a:b);
%     count = count + b-a+1;
% end
% h1 = histogram(temp1,100,'Normalization','pdf','DisplayStyle','stairs');
% N1 = h1.Values;
% Edges = h1.BinEdges;
% hold on;

% load('C:\Users\Administrator\Desktop\result_total_nondeg_p4_h2_res100_num100.mat')
% count = 1;
% for i = 1:num
%     len = length(q_store{4,i});
%     a = round(0.45*len);
%     b = round(0.55*len);
%     temp2(count:count+b-a) = q_store{4,i}(a:b);
%     count = count + b-a+1;
% end
% h2 = histogram(temp2,100,'Normalization','pdf','DisplayStyle','stairs');
% N2 = h2.Values;
% legend('h = 0.1','h = 2')
% xlabel('r')
%ylabel('count*binwidth')

% axis([0.45 0.55 0.06 0.24])
% axis([0.45 0.55 0.06 0.35])

% figure;
load('C:\Users\Administrator\Desktop\main_N26_p4_h0.1_No1.mat')
x = (e-e(1))/(e(end)-e(1));
for i = 1:NN
    if x(i) >0.45
        a = i;
        break
    end
end

for i = a:NN
    if x(i) >0.55
        b = i;
        break
    end
end
temp = zeros(b-a+1,4);
temp(1:b-a+1,1) = x(a:b);
temp(1:b-a+1,2) = mm(1,a:b)';
% plot(x,mm(1,:))
% hold on
load('C:\Users\Administrator\Desktop\main_N22_p4_h2_No1.mat')
x = (e-e(1))/(e(end)-e(1));
for i = 1:NN
    if x(i) >0.45
        a = i;
        break
    end
end

for i = a:NN
    if x(i) >0.55
        b = i;
        break
    end
end
temp(1:b-a+1,3) = x(a:b);
temp(1:b-a+1,4) = mm(1,a:b)';
% temp(b-a+2:end,3) = [];
% temp(b-a+2:end,4) = [];
% plot(x,mm(1,:))
% axis([0.45 0.55 0.06 0.35])

save('p=4,h=0_1,order_22,26.txt','temp','-ascii');



