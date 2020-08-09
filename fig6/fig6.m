clc, clear all
clf;
name = 'fig6_simulation_result_100_replica_with_TEs';
figname = sprintf('x%s',name);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');

clear tmp1;
clear tmp;
clear lys;
clear inf;
clear ben;
clear neu;
clear proAve;
cd genemeansA
for i = 1:100
BaseName='genemeans_A';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    lys(i,:) = tmp(2,:);
    inf(i,:) = tmp(3,:);
    neu(i,:) = tmp(4,:);
    ben(i,:) = tmp(5,:);
    subplot(2,3,1)
    plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 3)
    hold on;
    plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 3)
    hold on;
    plot(tmp(1,:), neu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 3)
    hold on;
    plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 3)
    hold on;
    fclose(fid);
end
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben);

a = tmp(1,:);

subplot(2,3,1)
p1 = plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 3)
    hold on;
p2 = plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 3)
hold on;
p3 = plot(tmp(1,:), neu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 3)
hold on;
p4 = plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 3)
hold on;
p = plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,3), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
legend([p1, p2, p3, p4, p],  {'Excision genes', 'Re-infection genes','Neutral genes', 'Beneficial genes', 'Mean'},'fontsize',16,'interpreter','latex');
legend('boxoff')
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.1, 0.99999, 0.0, 0.0], 'String', "A", 'fontweight','bold', 'FontSize', 18)
cd ..

cd ismeansD
for i = 1:100
BaseName='ismeans_D';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp1 = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    islys(i,:) = tmp1(2,:);
    isinf(i,:) = tmp1(3,:);
    isneu(i,:) = tmp1(4,:);
    isben(i,:) = tmp1(5,:);
    b = mod(i, 5);
    if b ==0
    subplot(2,3,4)
    plot(tmp1(1,:), islys(i, :), 'color', 'r', 'LineWidth', 0.2)
    hold on;
    plot(tmp1(1,:), isinf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    plot(tmp1(1,:), isneu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 0.2)
    hold on;
    plot(tmp1(1,:), isben(i, :), 'color', 'g', 'LineWidth', 0.2)
    hold on;
    end
    fclose(fid);
end
proAve(:,1) = mean(islys);
proAve(:,2) = mean(isinf);
proAve(:,3) = mean(isneu);
proAve(:,4) = mean(isben);

a = tmp(1,:);
subplot(2,3,4)
    p1 = plot(tmp1(1,:), islys(i, :), 'color', 'r', 'LineWidth', 0.2)
    hold on;
    p2 = plot(tmp1(1,:), isinf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    p3 = plot(tmp1(1,:), isneu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 0.2)
    hold on;
    p4 = plot(tmp1(1,:), isben(i, :), 'color', 'g', 'LineWidth', 0.2)
    hold on;
p = plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,3), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
hold on;
legend([p1, p2, p3, p4, p],  {'Disruptions in excision genes', 'Disruptions in re-infection genes','Disruptions in neutral genes', 'Disruptions in beneficial genes', 'Mean'},'fontsize',16,'interpreter','latex');
legend('boxoff')
ylim([0 1.3]);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.1, 0.5, 0.0, 0.0], 'String', "D", 'fontweight','bold', 'FontSize', 18)
cd ..


clear tmp1;
clear tmp;
clear lys;
clear inf;
clear ben;
clear neu;
clear proAve;
cd genemeansB
for i = 1:100
BaseName='genemeans_B';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    lys(i,:) = tmp(2,:);
    inf(i,:) = tmp(3,:);
    neu(i,:) = tmp(4,:);
    ben(i,:) = tmp(5,:);
    subplot(2,3,2)
    plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 3)
    hold on;
    plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 3)
    hold on;
    plot(tmp(1,:), neu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 3)
    hold on;
    plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 3)
    hold on;
    fclose(fid);
end
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben);

a = tmp(1,:);

subplot(2,3,2)
p1 = plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 3)
    hold on;
p2 = plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 3)
hold on;
p3 = plot(tmp(1,:), neu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 3)
hold on;
p4 = plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 3)
hold on;
p = plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,3), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.38, 0.99999, 0.0, 0.0], 'String', "B", 'fontweight','bold', 'FontSize', 18)
cd ..

cd ismeansE
for i = 1:100
BaseName='ismeans_E';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp1 = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    islys(i,:) = tmp1(2,:);
    isinf(i,:) = tmp1(3,:);
    isneu(i,:) = tmp1(4,:);
    isben(i,:) = tmp1(5,:);
    b = mod(i, 5);
    if b ==0
    subplot(2,3,5)
    plot(tmp1(1,:), islys(i, :), 'color', 'r', 'LineWidth', 0.2)
    hold on;
    plot(tmp1(1,:), isinf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    plot(tmp1(1,:), isneu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 0.2)
    hold on;
    plot(tmp1(1,:), isben(i, :), 'color', 'g', 'LineWidth', 0.2)
    hold on;
    end
    fclose(fid);
end
proAve(:,1) = mean(islys);
proAve(:,2) = mean(isinf);
proAve(:,3) = mean(isneu);
proAve(:,4) = mean(isben);

a = tmp(1,:);
subplot(2,3,5)
    p1 = plot(tmp1(1,:), islys(i, :), 'color', 'r', 'LineWidth', 0.2)
    hold on;
    p2 = plot(tmp1(1,:), isinf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    p3 = plot(tmp1(1,:), isneu(i, :), 'color', [0.75, 0.75, 0.75], 'LineWidth', 0.2)
    hold on;
    p4 = plot(tmp1(1,:), isben(i, :), 'color', 'g', 'LineWidth', 0.2)
    hold on;
p = plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,3), 'k', 'LineWidth', 1);
hold on;

plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
hold on;
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.38, 0.5, 0.0, 0.0], 'String', "E", 'fontweight','bold', 'FontSize', 18)
cd ..

clear tmp1;
clear tmp;
clear lys;
clear inf;
clear ben;
clear neu;
clear proAve;
    tmp1 = zeros(5, 1000);
    lys = zeros(100, 1000);
    inf = zeros(100, 1000);
    neu = zeros(100, 1000);
    ben = zeros(100, 1000);
% a = linspace(0, 200, 1000);
last = 400;
cd genemeansC
for i=1:100
    BaseName='genemeans_C';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    a = tmp(1,:);
    tmp1(:, 1:length(tmp)) = tmp;
    lys(i,1:length(tmp)) = tmp1(2,1:length(tmp));
    inf(i,1:length(tmp)) = tmp1(3,1:length(tmp));
    neu(i, 1:length(tmp)) = tmp1(4,1:length(tmp));
    ben(i,1:length(tmp)) = tmp1(5,1:length(tmp));
     b = mod(i, 5);
    subplot(2,3,3)
    plot(a(1:last), lys(i, 1:last), 'color', 'r', 'LineWidth', 3)
    hold on;
    plot(a(1:last), inf(i, 1:last), 'color', 'y', 'LineWidth', 3)
    hold on;
     plot(a(1:last), neu(i, 1:last), 'color', [0.75, 0.75, 0.75], 'LineWidth', 3)
    hold on;
    plot(a(1:last), ben(i, 1:last), 'color', 'g', 'LineWidth', 3)
    hold on;
    fclose(fid);
end
    
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben);

a = tmp(1,:);

subplot(2,3,3)
plot(a(1:last), proAve(1:last,1), 'k', 'LineWidth', 1);
hold on;

plot(a(1:last), proAve(1:last,2), 'k','LineWidth', 1);
hold on;

plot(a(1:last), proAve(1:last,3), 'k', 'LineWidth', 1);
hold on;

plot(a(1:last), proAve(1:last,4), 'k', 'LineWidth', 1);
hold on
ylim([0 8]);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.67, 0.9999, 0.0, 0.0], 'String', "C", 'fontweight','bold', 'FontSize', 18)
cd ..

clear istmp1;
clear istmp;
clear islys;
clear isinf;
clear isben;
clear isneu;
clear proAve;
istmp1 = zeros(5, 1000);
islys = zeros(100, 1000);
isinf = zeros(100, 1000);
isneu = zeros(100, 1000);
isben = zeros(100, 1000);
cd ismeansF
for i=1:100
    BaseName='ismeans_F';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    istmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    a = istmp(1,:);
    istmp1(:, 1:length(istmp)) = istmp;
    islys(i,1:length(istmp)) = istmp1(2,1:length(istmp));
    isinf(i,1:length(istmp)) = istmp1(3,1:length(istmp));
    isneu(i, 1:length(istmp)) = istmp1(4,1:length(istmp));
    isben(i,1:length(istmp)) = istmp1(5,1:length(istmp));
    c = mod(i,5);
    if c == 0
    subplot(2,3,6)
    plot(a(1:last), islys(i, 1:last), 'color', 'r', 'LineWidth', 3)
    hold on;
    plot(a(1:last), isinf(i, 1:last), 'color', 'y', 'LineWidth', 3)
    hold on;
    plot(a(1:last), isneu(i, 1:last), 'color', [0.75, 0.75, 0.75], 'LineWidth', 3)
    hold on;
    plot(a(1:last), isben(i, 1:last), 'color', 'g', 'LineWidth', 3)
    hold on;
    end
    fclose(fid);
end
proAve(:,1) = mean(islys);
proAve(:,2) = mean(isinf);
proAve(:,3) = mean(isneu);
proAve(:,4) = mean(isben);

a = tmp(1,:);
subplot(2,3,6)
 plot(a(1:last), proAve(1:last,1), 'k', 'LineWidth', 1);
hold on;

plot(a(1:last), proAve(1:last,2), 'k','LineWidth', 1);
hold on;
plot(a(1:last), proAve(1:last,3), 'k', 'LineWidth', 1);
hold on;

plot(a(1:last), proAve(1:last,4), 'k', 'LineWidth', 1);

xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.67, 0.5, 0.0, 0.0], 'String', "F", 'fontweight','bold', 'FontSize', 18)
cd ..

stampname(0.959468176914779, 0.0759625390218522,90);
creat(name);

