clc, clear all
clf;
name = 'fig8_trans';
figname = sprintf('x%s',name);
set(gcf,'position', [647 505 1202 457]);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');


subplot(1,2,1)
load totalsbylength_data1.txt
t = totalsbylength_data1;   % row 1 is non-IS, row 2 is IS transposases
ls = 5:5:70;                % length bins

tnorm(1,:) = t(1,:)./ls;    % normalize to average number of TEs per kB
tnorm(2,:) = t(2,:)./ls;
H = bar(ls,tnorm');
H(1).FaceColor = [0.2 0.2 0.2];
H(2).FaceColor = [0.9 0.9 0.9];
annotation('textbox', [0.06, 0.999999, 0.0, 0.0], 'String', "A", 'fontweight','bold', 'FontSize', 18)

xlabel('Prophage length (kbp)', 'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel('Mean transposes per kbp', 'fontsize',18,'verticalalignment','bottom','interpreter','latex');
axis([2.5 72.5 0 4])

%legend('non-IS transposases','IS transposases','location','NE');
legend('Non-IS transposases','IS transposases','location','NE', 'interpreter','latex', 'fontsize',16);
set(gca,'fontsize',16)
set(gcf,'color','w');

subplot(1,2,2)
load totalsbylength_data2.txt
t = totalsbylength_data2;
ls = 5:5:70;

tnorm(1,:) = t(1,:)./ls;
tnorm(2,:) = t(2,:)./ls;

H = bar(ls,tnorm');
H(1).FaceColor = [0.2 0.2 0.2];
H(2).FaceColor = [0.9 0.9 0.9];

xlabel('Prophage length (kbp)', 'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel('Mean transposes per kbp', 'fontsize',18,'verticalalignment','bottom','interpreter','latex','fontsize',16);
axis([2.5 72.5 0 8])
annotation('textbox', [0.5, 0.999999, 0.0, 0.0], 'String', "B", 'fontweight','bold', 'FontSize', 18)
set(gca,'fontsize',16)

stampname(81.0422662146837, -0.620681563524794,90);
creat(name);

