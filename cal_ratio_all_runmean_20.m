diro='./CAL/';
R='NCEP';
SLA='EGA';
syear=70;

filo1=[ diro 'cal_RCP26_' R '_' SLA '_runmean_20.nc'];
filo2=[ diro 'cal_RCP45_' R '_' SLA '_runmean_20.nc'];
filo3=[ diro 'cal_RCP60_' R '_' SLA '_runmean_20.nc'];
filo4=[ diro 'cal_RCP85_' R '_' SLA '_runmean_20.nc'];

cal1=ncread(filo1,'cal');
cal2=ncread(filo2,'cal');
cal3=ncread(filo3,'cal');
cal4=ncread(filo4,'cal');

nmodel1=length(ncread(filo1,'model'));
nmodel2=length(ncread(filo2,'model'));
nmodel3=length(ncread(filo3,'model'));
nmodel4=length(ncread(filo4,'model'));

f1=ncread(filo1,'f')';
f2=ncread(filo2,'f')';
f3=ncread(filo3,'f')';
f4=ncread(filo4,'f')';

%-----------------------------------------------------------------------------
% calculate mean and STD
%-----------------------------------------------------------------------------

byear=2020;
eyear=2099;

byear=byear-1948-syear+1;
eyear=eyear-1948-syear+1;

nmean=eyear-byear+1; %eyear-byear+1;

w1=repmat(ncread(filo1,'weight'),[1 nmean])';
w2=repmat(ncread(filo2,'weight'),[1 nmean])';
w3=repmat(ncread(filo3,'weight'),[1 nmean])';
w4=repmat(ncread(filo4,'weight'),[1 nmean])';

p1=sum(w1.*f1((syear+byear):(syear+eyear),:),2);
p2=sum(w2.*f2((syear+byear):(syear+eyear),:),2);
p3=sum(w3.*f3((syear+byear):(syear+eyear),:),2);
p4=sum(w4.*f4((syear+byear):(syear+eyear),:),2);
  
sigma1=mean(sqrt(sum(w1.*(f1((syear+byear):(syear+eyear),:)-repmat(p1,[1 nmodel1])).^2,2)));
sigma2=mean(sqrt(sum(w2.*(f2((syear+byear):(syear+eyear),:)-repmat(p2,[1 nmodel2])).^2,2)));
sigma3=mean(sqrt(sum(w3.*(f3((syear+byear):(syear+eyear),:)-repmat(p3,[1 nmodel3])).^2,2)));
sigma4=mean(sqrt(sum(w4.*(f4((syear+byear):(syear+eyear),:)-repmat(p4,[1 nmodel4])).^2,2)));
  
p1_avg=mean(f1((syear+byear):(syear+eyear),:),2);
p2_avg=mean(f2((syear+byear):(syear+eyear),:),2);
p3_avg=mean(f3((syear+byear):(syear+eyear),:),2);
p4_avg=mean(f4((syear+byear):(syear+eyear),:),2);

sigma1_avg=mean(sqrt(mean((f1((syear+byear):(syear+eyear),:)-repmat(p1_avg,[1 nmodel1])).^2,2)));
sigma2_avg=mean(sqrt(mean((f2((syear+byear):(syear+eyear),:)-repmat(p2_avg,[1 nmodel2])).^2,2)));
sigma3_avg=mean(sqrt(mean((f3((syear+byear):(syear+eyear),:)-repmat(p3_avg,[1 nmodel3])).^2,2)));
sigma4_avg=mean(sqrt(mean((f4((syear+byear):(syear+eyear),:)-repmat(p4_avg,[1 nmodel4])).^2,2)));
  
c1=(((squeeze(cal1(1,1,2:end))+squeeze(cal1(1,2,2:end)))*sigma1)./((2*icdf('norm',1-squeeze(sum(cal1(3,:,2:end)))/2,0,1))*sigma1_avg));
c2=(((squeeze(cal2(1,1,2:end))+squeeze(cal2(1,2,2:end)))*sigma2)./((2*icdf('norm',1-squeeze(sum(cal2(3,:,2:end)))/2,0,1))*sigma2_avg));
c3=(((squeeze(cal3(1,1,2:end))+squeeze(cal3(1,2,2:end)))*sigma3)./((2*icdf('norm',1-squeeze(sum(cal3(3,:,2:end)))/2,0,1))*sigma3_avg));
c4=(((squeeze(cal4(1,1,2:end))+squeeze(cal4(1,2,2:end)))*sigma4)./((2*icdf('norm',1-squeeze(sum(cal4(3,:,2:end)))/2,0,1))*sigma4_avg));

S1=repmat(sigma1./sigma1_avg,9,1);
S2=repmat(sigma2./sigma2_avg,9,1);
S3=repmat(sigma3./sigma3_avg,9,1);
S4=repmat(sigma4./sigma4_avg,9,1);

g1=(((squeeze(cal1(1,1,2:end))+squeeze(cal1(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal1(3,:,2:end)))/2,0,1))));
g2=(((squeeze(cal2(1,1,2:end))+squeeze(cal2(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal2(3,:,2:end)))/2,0,1))));
g3=(((squeeze(cal3(1,1,2:end))+squeeze(cal3(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal3(3,:,2:end)))/2,0,1))));
g4=(((squeeze(cal4(1,1,2:end))+squeeze(cal4(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal4(3,:,2:end)))/2,0,1))));


x1=1-(squeeze(sum(cal1(3,:,2:end))));
x2=1-(squeeze(sum(cal2(3,:,2:end))));
x3=1-(squeeze(sum(cal3(3,:,2:end))));
x4=1-(squeeze(sum(cal4(3,:,2:end))));

YL=2;
YD=0.2;
YS=0;

r=0.75;

h3=figure(3);
set(gcf,'Resize','off')
set(h3,'position',[100 100 1000 600]);
dx=0.43;
dy=0.38;
ddx=(1.00-2*dx)/3;
ddy=(0.92-2*dy)/2;

FontSize=16;
FontSizel=16;

s1=subplot('position',[0*dx+1*ddx+0.01 1*dy+2.5*ddy dx dy]);
i1=plot(x1',[g1' ; S1' ; c1' ; repmat(1,1,9)],'LineWidth',2);
set(i1(4),'Color',[0 0 0])
set(gca,'FontSize',FontSize)
title('RCP2.6')
xlim([0 1]);
ylim([YS YL]);
ylabel('Ratio')
%xlabel('Confidence level')
l=legend('(\gamma_u+\gamma_d)/(2\delta_G)','STD ratio','Interval ratio','Unity line');
set(l,'FontSize',FontSizel)
set(l,'Location','NorthWest')
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

s2=subplot('position',[1*dx+2*ddx+0.02 1*dy+2.5*ddy dx dy]);
i2=plot(x2',[g2' ; S2' ; c2' ; repmat(1,1,9)],'LineWidth',2);
set(i2(4),'Color',[0 0 0])
set(gca,'FontSize',FontSize)
title('RCP4.5')
xlim([0 1]);
ylim([YS YL]);
%ylabel('Ratio')
%xlabel('Confidence level')
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

s3=subplot('position',[0*dx+1*ddx+0.01 0*dy+1.*ddy dx dy]);
i3=plot(x3',[g3' ; S3' ; c3' ; repmat(1,1,9)],'LineWidth',2);
set(i3(4),'Color',[0 0 0])
set(gca,'FontSize',FontSize)
title('RCP6.0')
xlim([0 1]);
ylim([YS YL]);
ylabel('Ratio')
xlabel('Confidence level')
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

s4=subplot('position',[1*dx+2*ddx+0.02 0*dy+1.*ddy dx dy]);
i4=plot(x4',[g4' ; S4' ; c4' ; repmat(1,1,9)],'LineWidth',2);
set(i4(4),'Color',[0 0 0])
set(gca,'FontSize',FontSize)
title('RCP8.5')
xlim([0 1]);
ylim([YS YL]);
%ylabel('Ratio')
xlabel('Confidence level')
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

fname=['PLOTS/cal_ratio_all_runmean_20'];
%print(h3,'-deps','-r300',fname);
savefig(h3,fname)
saveas(h3,fname,'eps2c')
