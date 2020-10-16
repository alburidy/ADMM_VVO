% Programmed by A. Alburidy and L. Fan
% arburidy@gmail.com
% If you find this code useful for your research, please cite our paper at:
% https://github.com/alburidy/ADMM-VVO-Optimization
%==========================
% Convergence prove graph
K = length(history.obj);
m = {'+';'o';'*';'.';'x';'s';'d';'^';'v';'>';'<';'p';'h';'o';'*';};
for k=1:SD.tpn
    legend_lambda{k,1}=sprintf('LTC lambda %d',k); %#ok<*SAGROW>
end
i1=0;
for k=SD.tpn+1:SD.tpn+SD.cbn
    i1=i1+1;
    legend_lambda{k,1}=sprintf('CB lambda %d',i1);
end

% figure(2)
% plot(1:K, history.obj, 'k', 'MarkerSize', 10, 'LineWidth', 2);
% title('The Problem''s Convergence Check');
% ylabel('Objective Function'); xlabel('iter (k)');
% grid on
% xlim([1 K])


figure(3)
theta_deg=rad2deg(theta);
theta1 = linspace(0,2*pi);
rho1 = repmat(0.9,length(theta1));
theta2 = linspace(0,2*pi);
rho2 = repmat(1.1,length(theta1));
polarplot(theta_deg,v,'+k', 'LineWidth', 1.5)
hold on
polarplot(theta1,rho1,'-r', 'LineWidth', 1.5)
hold on
polarplot(theta2,rho2,'-r', 'LineWidth', 1.5)
hold off
title('Nodal Voltage Magnetude and Angles');

% figure(4)
% subplot(3,1,1);
% plot(1:K, history.obj, 'k', 'MarkerSize', 10, 'LineWidth', 2);
% title('Sub-Problem 1 Objective Function''s Value');
% ylabel('Objective Function'); xlabel('iter (k)');
% grid on
% xlim([1 K])
% subplot(3,1,2);
% a=plot(1:K, history.test_1, 'k',1:K, repmat(epison_1,K,1), 'r--',  'LineWidth', 1);
% grid on
% title('Primal Residual Convergence Check');
% ylabel('||r||_2');
% ylim([min(history.test_1) max(history.test_1)])
% xlim([1 K])
% legend('||r||_2','\epsilon^{pri}')
% subplot(3,1,3);
% b=plot(1:K, history.test_2, 'k',1:K,repmat(epison_2,K,1), 'r--', 'LineWidth', 1);
% grid on
% title('Dual Residual Convergence Check');
% ylabel('||s||_2'); xlabel('iter (k)');
% ylim([min(history.test_2-epison_2)*0.9 max(history.test_2+epison_2)*1.1])
% xlim([1 K])
% legend('||s||_2','\epsilon^{dual}')


% figure(5)
% subplot(2,1,1);
% AA = plot(1:K, history.lambda,'LineWidth', 1.5);
% set(AA,{'Marker'},m(1:size(history.lambda,2)));
% set(AA, {'color'}, num2cell(colors, 2));
% grid on
% title('Evolution of Scaled Lagrangian Multiplier');
% ylabel('\lambda');
% xlim([1 K])
% legend(legend_lambda)
% subplot(2,1,2);
% PLOT = plot(1:itr,history.convergence,'LineWidth', 1.5);
% set(PLOT, {'color'}, num2cell(colors, 2));
% grid on;
% ylabel('h(x,u)');
% legend(legends);
% title('Convergence Progress'); xlabel('iteration #');

% figure(6)
% subplot(2,2,1);
% a=plot(1:K, history.test_1, 'k',1:K, repmat(epison_1,K,1), 'r--',  'LineWidth', 1);
% grid on
% title('(a)', 'FontSize', 12);
% xlabel('iteration (k)');
% xlim([1 K])
% ylim([-0.02 0.2])
% legend('||r||_2','\epsilon^{pri}')
% subplot(2,2,2);
% b=plot(1:K, history.test_2, 'k',1:K,repmat(epison_2,K,1), 'r--', 'LineWidth', 1);
% grid on
% title('(b)', 'FontSize', 12);
% xlabel('iteration (k)');
% xlim([1 K])
% ylim([-100 1000])
% legend('||s||_2','\epsilon^{dual}')
% subplot(2,2,3);
% plot(1:K, history.obj, 'k', 'MarkerSize', 10, 'LineWidth', 2);
% title('(c)', 'FontSize', 12);
% legend('f_x(x)'); xlabel('iteration (k)');
% grid on
% xlim([1 K])
% subplot(2,2,4);
% plot(1:K, history.rho, 'k', 'MarkerSize', 10, 'LineWidth', 2);
% title('(d)', 'FontSize', 12);
% xlabel('iteration (k)');
% grid on
% legend('\rho')
% xlim([1 K])
% ylim([0 250])


figure(7)
subplot(1,2,1);
AA = plot(1:K, history.lambda,'LineWidth', 1.5);
set(AA,{'Marker'},m(1:size(history.lambda,2)));
grid on
title('(a)', 'FontSize', 10);
ylabel('\lambda');
xlim([1 K])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.3f'))
xlabel('iteration (k)');
legend('$LTC_{1}$','$LTC_{2}$','$LTC_{3}$','$SCB_{1}$','$SCB_{2}$','$SCB_{3}$','$SCB_{4}$', 'interpreter', 'latex');
% ylim([-0.006 0.1])
subplot(1,2,2);
PLOT = plot(1:itr,history.convergence,'LineWidth', 1.5);
set(PLOT,{'Marker'},m(1:size(history.lambda,2)));
grid on;
xlim([1 K])
ylabel('h(x,u)');
title('(b)', 'FontSize', 10);
xlabel('iteration (k)');
% ylim([-0.13 0.02])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

% ========================================
figure(6)
subplot(2,2,1);
a=plot(1:K, history.test_1, 'k',1:K, repmat(epison_1,K,1), 'r--',  'LineWidth', 1.5);
grid on
title('(a) Primal Residual', 'FontSize', 10);
xlabel('iteration (k)');
xlim([1 K])
ylim([-0.02 0.15])
legend('||r||_2','\epsilon^{pri}')
% Specify x-label position
Xlb=mean(xlim);
Ylb=-0.038;
xlabel('iteration (k)','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 
subplot(2,2,2);
b=plot(1:K, history.test_2, 'k',1:K,repmat(epison_2,K,1), 'r--', 'LineWidth', 1.5);
grid on
title('(b) Dual Residual', 'FontSize', 10);
xlabel('iteration (k)');
xlim([1 K])
ylim([-50 350])
legend('||s||_2','\epsilon^{dual}');
% Specify x-label position
Ylb=-90;
xlabel('iteration (k)','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center');
subplot(2,2,3);
plot(1:K, history.obj, 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
title('(c) Objective Function Value', 'FontSize', 10);
legend('f_x(x)'); xlabel('iteration (k)');
grid on
xlim([1 K])
% ylim([-5 30])
% Specify x-label position
Ylm=ylim;
Ylb=-0.037;
xlabel('iteration (k)','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center');
subplot(2,2,4);
plot(1:K, history.rho, 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
title('(d) Penalty Parameter Value', 'FontSize', 10);
grid on
legend('\rho')
xlim([1 K]);
Xlb=mean(xlim);
ylim([0 250]);
% Specify x-label position
Ylm=ylim;
Ylb=-20;
xlabel('iteration (k)','Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center'); 