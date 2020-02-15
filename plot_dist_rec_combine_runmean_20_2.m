load CAL/calc_pdf_rand_ALL_2099-2099_runmean_gaus_20.mat

PM1_05g=PM1_05;
PM2_05g=PM2_05;
PM3_05g=PM3_05;
PM4_05g=PM4_05;
sg=s;
load CAL/calc_pdf_rand_ALL_2099-2099_runmean_20.mat
XMaa=XM;
PM1a=PM1;
PM2a=PM2;
PM3a=PM3;
PM4a=PM4;
X1a=x1(1:end-1,1)+diff(x1(:,1))/2;
X2a=x2(1:end-1,1)+diff(x2(:,1))/2;
X3a=x3(1:end-1,1)+diff(x3(:,1))/2;
X4a=x4(1:end-1,1)+diff(x4(:,1))/2;
Y1a=PI1;Y2a=PI2;Y3a=PI3;Y4a=PI4;

Q2=Q;
S2=S2
FS=12;
s1=[PM1_05 PM1_05g 1.0 PM2_05 PM2_05g 1.8 PM3_05 PM3_05g 2.2 PM4_05 PM4_05g 3.7 ];
ii1=[1 4 7 10];
ii2=[2 5 8 11];
ii3=[3 6 9 12];
s2(:,ii1)=s;
s2(:,ii2)=sg;
s2(1,ii3)=[0.3 1.1 1.4 2.6];
s2(2,ii3)=[0.6 1.3 1.7 3.0];
s2(3,ii3)=[1.4 2.3 2.7 4.4];
s2(4,ii3)=[1.7 2.6 3.1 4.8];

s2=s2([2 3 1 4],:);

data = rand(100,12);
h    = boxplot(data,'BoxStyle','filled','MedianStyle','target','Widths',1);
set(gca,'ylim',[0 5])
set(gca,'FontSize',FS)
set(gca,'XTick',2:3:11)
set(gca,'XTickLabel',{'RCP2.6','RCP4.5','RCP6.0','RCP8.5'})
ylabel('GMST change [^{\circ}C]')
set(h(1,:),{'Ydata'},num2cell(s2(3:4,:),1)')
set(h(2,:),{'Ydata'},num2cell(s2(2:-1:1,:),1)')
set(h(3,:),{'Ydata'},num2cell(s1,1)')
set(h(4,:),{'Ydata'},num2cell(s1,1)')
set(h(:,ii2),'Color',[0.8500, 0.3250, 0.0980])
set(h(:,ii1),'Color',[0.4660, 0.6740, 0.1880])
set(h(:,ii3),'Color',[0, 0.4470, 0.7410])
set(h(:,ii2),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
set(h(:,ii1),'MarkerEdgeColor',[0.4660, 0.6740, 0.1880])
set(h(:,ii3),'MarkerEdgeColor',[0, 0.4470, 0.7410])
l=legend(h(1,1:3),'AR+EGA','GAUS+AVG','IPCC');
set(l,'Location','NorthWest');
set(h(1,:),'LineWidth',2)
set(h(2,:),'LineWidth',6)
title('2080-2099')

fname=['PLOTS/box_year_2080-2099_runmean_20'];
print(gcf,'-dpng','-r300',fname);
savefig(gcf,fname)

load CAL/calc_pdf_rand_ALL_2065-2065_runmean_gaus_20.mat
PM1_05g=PM1_05;
PM2_05g=PM2_05;
PM3_05g=PM3_05;
PM4_05g=PM4_05;
sg=s;
load CAL/calc_pdf_rand_ALL_2065-2065_runmean_20.mat
XMbb=XM;
PM1b=PM1;
PM2b=PM2;
PM3b=PM3;
PM4b=PM4;

X1b=x1(1:end-1,1)+diff(x1(:,1))/2;
X2b=x2(1:end-1,1)+diff(x2(:,1))/2;
X3b=x3(1:end-1,1)+diff(x3(:,1))/2;
X4b=x4(1:end-1,1)+diff(x4(:,1))/2;
Y1b=PI1;Y2b=PI1;Y3b=PI3;Y4b=PI4;

Q1=Q;
S1=S2
s1=[PM1_05 PM1_05g 1.0 PM2_05 PM2_05g 1.4 PM3_05 PM3_05g 1.3 PM4_05 PM4_05g 2.0 ];
ii1=[1 4 7 10];
ii2=[2 5 8 11];
ii3=[3 6 9 12];
s2(:,ii1)=s;
s2(:,ii2)=sg;
s2(1,ii3)=[0.4 0.9 0.8 1.4];
s2(2,ii3)=[0.7 1.1 1.0 1.6];
s2(3,ii3)=[1.3 1.7 1.6 2.4];
s2(4,ii3)=[1.6 2.0 1.8 2.6];

s2=s2([2 3 1 4],:);

data = rand(100,12);
h    = boxplot(data,'BoxStyle','filled','MedianStyle','target','Widths',1);
set(gca,'ylim',[0 5])
set(gca,'FontSize',FS)
set(gca,'XTick',2:3:11)
set(gca,'XTickLabel',{'RCP2.6','RCP4.5','RCP6.0','RCP8.5'})
ylabel('GMST change [^{\circ}C]')
set(h(1,:),{'Ydata'},num2cell(s2(3:4,:),1)')
set(h(2,:),{'Ydata'},num2cell(s2(2:-1:1,:),1)')
set(h(3,:),{'Ydata'},num2cell(s1,1)')
set(h(4,:),{'Ydata'},num2cell(s1,1)')
set(h(:,ii2),'Color',[0.8500, 0.3250, 0.0980])
set(h(:,ii1),'Color',[0.4660, 0.6740, 0.1880])
set(h(:,ii3),'Color',[0, 0.4470, 0.7410])
set(h(:,ii2),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
set(h(:,ii1),'MarkerEdgeColor',[0.4660, 0.6740, 0.1880])
set(h(:,ii3),'MarkerEdgeColor',[0, 0.4470, 0.7410])
l=legend(h(1,1:3),'AR+EGA','GAUS+AVG','IPCC');
set(l,'Location','NorthWest');
set(h(1,:),'LineWidth',2)
set(h(2,:),'LineWidth',6)
title('2046-2065')

fname=['PLOTS/box_year_2046-2065_runmean_20'];
print(gcf,'-dpng','-r300',fname);
savefig(gcf,fname)








