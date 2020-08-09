clc, clear all
clf;
name = 'fig4_simulations_100_replica';
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
tmp1 = zeros(5, 1000);
lys = zeros(100, 1000);
inf = zeros(100, 1000);
neu = zeros(100, 1000);
ben = zeros(100, 1000);
% a = linspace(0, 200, 1000);
last = 400;
cd genemeansA
for i=1:100
    BaseName='genemeans_A';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    a = tmp(1,:);
    tmp1(:, 1:length(tmp)) = tmp;
    lys(i,1:length(tmp)) = tmp1(2,1:length(tmp));
    inf(i,1:length(tmp)) = tmp1(3,1:length(tmp));
    neu(i, 1:length(tmp)) = tmp1(4,1:length(tmp));
    ben(i,1:length(tmp)) = tmp1(5,1:length(tmp));
    b = mod(i,5);
    if b == 0
    subplot(2,2,1)
    plot(a(1:last), lys(i, 1:last), 'color', 'r', 'LineWidth', 2)
    hold on;
    plot(a(1:last), inf(i, 1:last), 'color', 'y', 'LineWidth', 02)
    hold on;
    plot(a(1:last), ben(i, 1:last), 'color', 'g', 'LineWidth', 02)
    hold on;
    end
        fclose(fid);
end
cd ..
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben); 
subplot(2,2,1)
p1 = plot(a(1:last), lys(end, 1:last), 'color', 'r', 'LineWidth',02)
hold on;
p2 = plot(a(1:last), inf(end, 1:last), 'color', 'y', 'LineWidth', 02)
hold on;
p3 = plot(a(1:last), ben(end, 1:last), 'color', 'g', 'LineWidth', 02)
hold on;
p = plot(a(1:last), proAve(1:last,1), 'k', 'LineWidth', 1);
hold on;
plot(a(1:last), proAve(1:last,2), 'k','LineWidth', 1);
hold on;
plot(a(1:last), proAve(1:last,4), 'k', 'LineWidth', 1);
hold on;
set(gca,'fontsize',16);
legend([p1, p2, p3, p],  {'Excision genes','Re-infection genes', 'Beneficial genes', 'Mean'},'fontsize',16,'interpreter','latex');
legend('boxoff')
xlim([0,40]);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.05, 0.999999, 0.0, 0.0], 'String', "A", 'fontweight','bold', 'FontSize', 18)
hold off;

clear tmp1;
clear tmp;
clear lys;
clear inf;
clear ben;
clear neu;
clear proAve;
cd genemeansB
for i=1:100
    BaseName='genemeans_B';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    lys(i,:) = tmp(2,:);
    inf(i,:) = tmp(3,:);
    neu(i,:) = tmp(4,:);
    ben(i,:) = tmp(5,:);
    b = mod(i,5);
    if b == 0
    subplot(2,2,2)
    plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 02)
    hold on;
    plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 02)
    hold on;
    end
    fclose(fid);
end
cd ..
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben);
a = tmp(1,:);
subplot(2,2,2)
plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;
plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;
plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
hold on;
set(gca,'fontsize',16);
xlim([0 250]);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.5, 0.999999, 0.0, 0.0], 'String', "B", 'fontweight','bold', 'FontSize', 18)
hold off;


clear tmp;
clear lys;
clear inf;
clear ben;
clear neu;
clear proAve;
cd genemeansC
for i=1:100
    BaseName='genemeans_C';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    lys(i,:) = tmp(2,:);
    inf(i,:) = tmp(3,:);
    neu(i,:) = tmp(4,:);
    ben(i,:) = tmp(5,:);
    b = mod(i,5);
    if b == 0
       
    subplot(2,2,3)
    plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 0.2)
    hold on;
    plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 02)
    hold on;
    end
    fclose(fid);
end
cd ..
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben);
a = tmp(1,:);
subplot(2,2,3)
plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;
plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;
plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
hold on;
set(gca,'fontsize',16);
xlim([0 1500]);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.05, 0.5, 0.0, 0.0], 'String', "C", 'fontweight','bold', 'FontSize', 18)
hold off;
 
clear tmp;
clear lys;
clear inf;
clear ben;
clear neu;
clear proAve;
cd genemeansD
for i=1:100
    BaseName='genemeans_D';
    FileName =[BaseName,num2str(i)];
    fid = fopen(FileName,'r');
    tmp = fscanf(fid,'%f %f %f %f %f', [5,Inf]);
    lys(i,:) = tmp(2,:);
    inf(i,:) = tmp(3,:);
    neu(i,:) = tmp(4,:);
    ben(i,:) = tmp(5,:);
    b = mod(i,5);
    if b == 0
       
    subplot(2,2,4)
    plot(tmp(1,:), lys(i, :), 'color', 'r', 'LineWidth', 0.2)
    hold on;
    plot(tmp(1,:), inf(i, :), 'color', 'y', 'LineWidth', 0.2)
    hold on;
    plot(tmp(1,:), ben(i, :), 'color', 'g', 'LineWidth', 02)
    hold on;
    end
    fclose(fid);
end
cd ..
proAve(:,1) = mean(lys);
proAve(:,2) = mean(inf);
proAve(:,3) = mean(neu);
proAve(:,4) = mean(ben);
a = tmp(1,:);
subplot(2,2,4)
plot(a(1:end), proAve(1:end,1), 'k', 'LineWidth', 1);
hold on;
plot(a(1:end), proAve(1:end,2), 'k','LineWidth', 1);
hold on;
plot(a(1:end), proAve(1:end,4), 'k', 'LineWidth', 1);
hold on;
set(gca,'fontsize',16);
xlim([0 1500]);
xlabel({'Time'},'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Average number of genes'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.5, 0.5, 0.0, 0.0], 'String', "D", 'fontweight','bold', 'FontSize', 18)
hold off;

stampname(1691.18259144089, -0.947122217710215  ,90);

creat(name);

