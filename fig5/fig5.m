clc, clear all
clf;
name = 'fig5_simulations_without_TE';
figname = sprintf('x%s',name);
set(gcf,'position', [161 435 1502 402]);
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');

nLys = 1;      % Number of lysis genes
nInf = 8;     % Number of infectious genes
nNeu = 0;     % Number of Neutral genes
nBen = 1;     % Number of beneficial genes
nGenes = nLys+nInf+ nNeu + nBen; % Number of genes on an individual prophage
nType = 5;   % types of genes on a prophage
x = 0:nGenes; %
psumdataSim = zeros(2,nType);
sumLysIntact = 0; %sum of lysis genes in intact prophages
sumInfnIntact = 0; %sum of infn genes in intact prophages
%sumNeuIntact = 0; %sum of Neutral genes in intact prophages
sumBenIntact = 0; %sum of Beneficial genes in intact prophages
sumTranIntact = 0; %sum of Transposases genes in intact prophages
sumLysIncomp = 0; %sum of lysis genes in incomplete prophages
sumInfnIncomp = 0; %sum of infn genes in incomplete prophages
%sumNeuIncomp = 0; %sum of Neutral genes in incomplete prophages
sumBenIncomp = 0; %sum of Beneficial genes in incomplete prophages
sumTranIncomp = 0; %sum of Transposases genes in incomplete prophages
Intact = 0;      % Intact genes
Incomplete = 0;  % questionable genes 
tmp1 = load('prosim.out');
tmp = tmp1;
tmp2 = tmp1;
for i=1:length(tmp)
    iid = find(tmp(i,:)==-1);
    T(i) = length(iid);
    tmp(i, iid) = 0;
end
for i=1:length(tmp2)
    iid = find(tmp2(i,:)==-1);
    tmp2(i, iid) = 1;
end
 S = sum(tmp2')';
 L = tmp(:, 1:nLys);
 I = sum(tmp(:, nLys+1:nLys+nInf)')';
 B =tmp(:, nLys+nInf+ nNeu+1:nLys+nInf+ nNeu + nBen);
for j = 1:length(tmp)
        if (L(j) == nLys && I(j) == nInf)
         Intact = Intact + 1;
         sumLysIntact = sumLysIntact+L(j);
         sumInfnIntact = sumInfnIntact+I(j);
         sumBenIntact =  sumBenIntact+B(j);
         sumTranIntact = sumTranIntact+T(j);
        else
         Incomplete = Incomplete + 1;
         sumLysIncomp = sumLysIncomp+L(j);
         sumInfnIncomp = sumInfnIncomp+I(j);
         sumBenIncomp =  sumBenIncomp+B(j);
         sumTranIncomp = sumTranIncomp+T(j);
        end            
end
    sumdataSim = [sumLysIntact sumLysIncomp; sumInfnIntact sumInfnIncomp; sumBenIntact sumBenIncomp];
for j =1:2
    pdataSim(j,:) = sumdataSim(:, j)./sum(sumdataSim(:, j));
end

fprintf('Percentage of intact genes is :'); 
fprintf('%i\n', 100*(Intact/(Intact+Incomplete)));
fprintf('Percentage of incomplete phages is :'); 
fprintf('%i\n', 100*(Incomplete/(Intact+Incomplete)));
ikeep = find(max(pdataSim)>0.001);
[deltas,order] = sort((pdataSim(1,ikeep)-pdataSim(2,ikeep))./pdataSim(1,ikeep),'descend');
labels =[{ 'Excision'; 'Re-infection'; 'Beneficial'}];

subplot(1,3,2)
grey = 0.5*[1 1 1];
bar(100*(pdataSim(2,ikeep(order))-pdataSim(1,ikeep(order)))./pdataSim(1,ikeep(order)), 'FaceColor',grey,'EdgeColor',grey,'LineWidth',1);
xlim([0.5, 3.5]);
ylabel({'Percent change'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
set(gca,'xticklabel',labels(ikeep(order),:));
set(gca,'fontsize',16);
annotation('textbox', [0.36, 0.999999, 0.0, 0.0], 'String', "B", 'fontweight','bold', 'FontSize', 18)

subplot(1,3,1)
b = bar(pdataSim');
b(1).FaceColor = [0.2 0.2 0.2];
b(2).FaceColor = [0.8 0.8 0.8];
xlim([0.7, 3.3]);
legend({'Intact', 'Incomplete'},'fontsize',16,'interpreter','latex');
set(gca,'xticklabel',labels(ikeep(order),:))
set(gca,'fontsize',16);

ylabel({'Frequency'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.09, 0.999999, 0.0, 0.0], 'String', "A", 'fontweight','bold', 'FontSize', 18)

%-----------------------------------------------------------------------------------------------------
    Lys = zeros(length(x),1);
    Infn = zeros(length(x),1);
    Ben = zeros(length(x),1);
    Tran = zeros(length(x),1);
    
    lyspro = [S, L];
    infpro = [S, I];
    benpro = [S, B];
    tranpro = [S T'];
    [proSize, jj , kk] = unique(S);  % to find average number of genes against length in an individual tial
    % Frequency of each type of genes in an indiviual trial
    meanLys = accumarray(kk, (1:numel(kk))' , [], @(x) mean(lyspro(x,2)))./(proSize);
    meanInf = accumarray(kk, (1:numel(kk))' , [], @(x) mean(infpro(x,2)))./(proSize);
    meanBen = accumarray(kk, (1:numel(kk))' , [], @(x) mean(benpro(x,2)))./(proSize);
    meanTran = accumarray(kk,(1:numel(kk))', [], @(x) mean(tranpro(x,2)))./proSize;
    a = find(ismember(x, unique(S)));
    Lys(a) = meanLys;
    Infn(a) = meanInf;
    Ben(a) = meanBen;
    Tran(a) = meanTran;
    FinalLys = Lys;
    FinalInf = Infn;
    FinalBen = Ben;
    FinalTran = Tran;
    
[np,xp] = hist(S, 0:1:max(S));
np = np ./trapz(xp, np);

subplot(1,3,3)
bar( xp, 4*np, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [1 1 1] );
hold on;
xlim([0 nGenes+1]);
p1 = plot(x(2:end), FinalLys(2:end),'s-r', 'MarkerFaceColor', 'r', 'LineWidth', 2);
hold on;
p2 = plot(x(2:end), FinalInf(2:end), 's-y', 'MarkerFaceColor', 'y','LineWidth', 2);
hold on;
p3 = plot(x(2:end), FinalBen(2:end), 's-g','MarkerFaceColor', 'g', 'LineWidth', 2);
hold on;
%p4 = plot(x(2:end), FinalTran(2:end), 's-k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
set(gca,'fontsize',16);
legend([p1,p2,p3],{'Excision Genes', 'Re-infection Genes', 'Beneficial Genes'},'Location', 'Best', 'fontsize',16,'interpreter','latex');
legend('boxoff')
xlabel({'Prophage length (number of genes)'}, 'fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel({'Frequency of each type of gene'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
annotation('textbox', [0.65, 0.999999, 0.0, 0.0], 'String', "C", 'fontweight','bold', 'FontSize', 18)

stampname(13.3296894409938, -0.0894936708860758,90);

creat(name);
      
