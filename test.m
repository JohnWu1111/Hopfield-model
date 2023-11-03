% x = NN;
% % y = phi4_all(20,:);
% figure
% % semilogy(x,y)
% % x = 1./x;
% plot(x,mean(sus_all(45:55,:)).*NN)

% figure
% x = N_a;
% y = mean(corre_in_all(45:55,:));
% plot(x,y,'-*')
% xlabel('N')
% %axis([0 0.065 0 0.2])

% le = cell(1, lN);
% for i = 1:lN
%     le{i} = strcat('N', num2str(NN(i)));
% end
% plot(1/res:1/res:1, sus_all(:,4:6)');
% xlabel('\epsilon')
% ylabel('\chi')
% legend(le(4:6));

% x = 1/res:1/res:1;
% y1 = sus_all(:,4);
% y2 = sus_all2(:,4);
% plot(x,y1,x,y2)

% c = 0;
% clear qq
% for i = 1:num
%     temp = q_store{3,i};
%     len = length(temp);
%     a = round(0.45*len);
%     b = round(0.55*len);
% %     a = 1;
% %     b = len;
%     len1 = length(temp(a:b));
%     temp1 = sort(temp(a:b));
%     if temp1(1) == 0
%         continue
%     end
%     qq(c+1:c+len1) = temp(a:b);
%     c = c+len1;
% end
% 
% figure;
% histogram(qq,1000,'Normalization','pdf','DisplayStyle','stairs');
% %axis([0 1 0.5 1.5])
% mean(qq)

% le = cell(1, lN);
% 
% for i = 1:lN
%     le{i} = strcat('N', num2str(N_a(i)));
% end
% 
% y = zeros(res,3);
% y(:,1) = sus_all(:,1);
% y(:,2) = sus_all(:,3);
% y(:,3) = sus_all(:,5);
% 
% figure;
% plot(1/res:1/res:1, y');
% xlabel('\epsilon')
% ylabel('\chi')
% legend(le);

figure;
for i = 1:lN
    N = N_a(i);
    for j = 1:num
        NN = NN_store(i,j);
        Dq2_it = Dq2_store{i,j};
        Dq2_temp1 = Dq2_it*log(NN);
        Dq2_temp = mean(Dq2_temp1(floor(0.45*NN):floor(0.55*NN)));
        loglog(NN,exp(-Dq2_temp),'*')
        hold on;
    end
end
