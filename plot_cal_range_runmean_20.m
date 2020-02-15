R='NCEP';
SLA='EGA';

syear=70;

diro='./CAL/'

filo1=[ diro 'cal_RCP26_' R '_' SLA '_runmean_20.nc'];
filo2=[ diro 'cal_RCP45_' R '_' SLA '_runmean_20.nc'];
filo3=[ diro 'cal_RCP60_' R '_' SLA '_runmean_20.nc'];
filo4=[ diro 'cal_RCP85_' R '_' SLA '_runmean_20.nc'];

ncep=ncread(filo1,'y');

clm=mean(ncep(58));

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

nmean=20;

w1=repmat(ncread(filo1,'weight'),[1 nmean])';
w2=repmat(ncread(filo2,'weight'),[1 nmean])';
w3=repmat(ncread(filo3,'weight'),[1 nmean])';
w4=repmat(ncread(filo4,'weight'),[1 nmean])';


%-----------------------------------------------------------------------------
% calculate mean and STD
%-----------------------------------------------------------------------------

for i=1:2

  if i==1
    byear=2065;
    eyear=2065;
  else
    byear=2099;
    eyear=2099;
  end

  byear=byear-1948-syear+1;
  eyear=eyear-1948-syear+1;

  mtype=1;

  switch mtype
    case 1
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
  
    case 2
     p1=mean(f1((syear+byear):(syear+eyear),:),1);
     p2=mean(f2((syear+byear):(syear+eyear),:),1);
     p3=mean(f3((syear+byear):(syear+eyear),:),1);
     p4=mean(f4((syear+byear):(syear+eyear),:),1);
  
     sigma1=sqrt(sum(w1(1,:).*(p1-repmat(mean(p1),[1 length(p1)])).^2));
     sigma2=sqrt(sum(w2(1,:).*(p2-repmat(mean(p2),[1 length(p2)])).^2));
     sigma3=sqrt(sum(w3(1,:).*(p3-repmat(mean(p3),[1 length(p3)])).^2));
     sigma4=sqrt(sum(w4(1,:).*(p4-repmat(mean(p4),[1 length(p4)])).^2));
  
     p1_avg=mean(f1((syear+byear):(syear+eyear),:),1);
     p2_avg=mean(f2((syear+byear):(syear+eyear),:),1);
     p3_avg=mean(f3((syear+byear):(syear+eyear),:),1);
     p4_avg=mean(f4((syear+byear):(syear+eyear),:),1);
  
     sigma1_avg=sqrt(mean((p1_avg-repmat(mean(p1_avg),[1 length(p1_avg)])).^2));
     sigma2_avg=sqrt(mean((p2_avg-repmat(mean(p2_avg),[1 length(p2_avg)])).^2));
     sigma3_avg=sqrt(mean((p3_avg-repmat(mean(p3_avg),[1 length(p3_avg)])).^2));
     sigma4_avg=sqrt(mean((p4_avg-repmat(mean(p4_avg),[1 length(p4_avg)])).^2));

    case 3
     p1=sum(w1.*f1((syear+byear):(syear+eyear),:),2);
     p2=sum(w2.*f2((syear+byear):(syear+eyear),:),2);
     p3=sum(w3.*f3((syear+byear):(syear+eyear),:),2);
     p4=sum(w4.*f4((syear+byear):(syear+eyear),:),2);

     sigma1=sqrt(mean(sum(w1.*(f1((syear+byear):(syear+eyear),:)-repmat(p1,[1 nmodel1])).^2,2)));
     sigma2=sqrt(mean(sum(w2.*(f2((syear+byear):(syear+eyear),:)-repmat(p2,[1 nmodel2])).^2,2)));
     sigma3=sqrt(mean(sum(w3.*(f3((syear+byear):(syear+eyear),:)-repmat(p3,[1 nmodel3])).^2,2)));
     sigma4=sqrt(mean(sum(w4.*(f4((syear+byear):(syear+eyear),:)-repmat(p4,[1 nmodel4])).^2,2)));

     p1_avg=mean(f1((syear+byear):(syear+eyear),:),2);
     p2_avg=mean(f2((syear+byear):(syear+eyear),:),2);
     p3_avg=mean(f3((syear+byear):(syear+eyear),:),2);
     p4_avg=mean(f4((syear+byear):(syear+eyear),:),2);

     sigma1_avg=sqrt(mean(mean((f1((syear+byear):(syear+eyear),:)-repmat(p1_avg,[1 nmodel1])).^2,2)));
     sigma2_avg=sqrt(mean(mean((f2((syear+byear):(syear+eyear),:)-repmat(p2_avg,[1 nmodel2])).^2,2)));
     sigma3_avg=sqrt(mean(mean((f3((syear+byear):(syear+eyear),:)-repmat(p3_avg,[1 nmodel3])).^2,2)));
     sigma4_avg=sqrt(mean(mean((f4((syear+byear):(syear+eyear),:)-repmat(p4_avg,[1 nmodel4])).^2,2)));  
  end
  
  c1=(((squeeze(cal1(1,1,2:end))+squeeze(cal1(1,2,2:end)))*sigma1)./((2*icdf('norm',1-squeeze(sum(cal1(3,:,2:end)))/2,0,1))*sigma1_avg))
  c2=(((squeeze(cal2(1,1,2:end))+squeeze(cal2(1,2,2:end)))*sigma2)./((2*icdf('norm',1-squeeze(sum(cal2(3,:,2:end)))/2,0,1))*sigma2_avg))
  c3=(((squeeze(cal3(1,1,2:end))+squeeze(cal3(1,2,2:end)))*sigma3)./((2*icdf('norm',1-squeeze(sum(cal3(3,:,2:end)))/2,0,1))*sigma3_avg))
  c4=(((squeeze(cal4(1,1,2:end))+squeeze(cal4(1,2,2:end)))*sigma4)./((2*icdf('norm',1-squeeze(sum(cal4(3,:,2:end)))/2,0,1))*sigma4_avg))

  r1=(squeeze(cal1(1,1,2:end))+squeeze(cal1(1,2,2:end)))*sigma1;
  r2=(squeeze(cal2(1,1,2:end))+squeeze(cal2(1,2,2:end)))*sigma2;
  r3=(squeeze(cal3(1,1,2:end))+squeeze(cal3(1,2,2:end)))*sigma3;
  r4=(squeeze(cal4(1,1,2:end))+squeeze(cal4(1,2,2:end)))*sigma4;

  r1_avg=(2*icdf('norm',1-squeeze(sum(cal1(3,:,2:end)))/2,0,1)*sigma1_avg);
  r2_avg=(2*icdf('norm',1-squeeze(sum(cal2(3,:,2:end)))/2,0,1)*sigma2_avg);
  r3_avg=(2*icdf('norm',1-squeeze(sum(cal3(3,:,2:end)))/2,0,1)*sigma3_avg);
  r4_avg=(2*icdf('norm',1-squeeze(sum(cal4(3,:,2:end)))/2,0,1)*sigma4_avg);

  g1=squeeze(cal1(1,1:2,2:end))';
  g2=squeeze(cal2(1,1:2,2:end))';
  g3=squeeze(cal3(1,1:2,2:end))';
  g4=squeeze(cal4(1,1:2,2:end))';
 
  if i==1
    R1=[r1 r2 r3 r4]
    R1_avg=[r1_avg r2_avg r3_avg r4_avg]
    G1=[g1 g2 g3 g4];
  else
    R2=[r1 r2 r3 r4]
    R2_avg=[r1_avg r2_avg r3_avg r4_avg]
    G2=[g1 g2 g3 g4];
  end   
end

Ratio=[R1 R2]./[R1_avg R2_avg];

tex_table = latex(vpa(sym(round([(0.9:-0.1:0.1)' R1 R2],2))))
tex_table = latex(vpa(sym(round([(0.9:-0.1:0.1)' R1_avg R2_avg],2))))
tex_table = latex(vpa(sym(round([(0.9:-0.1:0.1)' Ratio],2))))
tex_table = latex(vpa(sym(round([(0.9:-0.1:0.1)' g1 g2 g3 g4],2))))

x1=1-(squeeze(sum(cal1(3,:,2:end))));
x2=1-(squeeze(sum(cal2(3,:,2:end))));
x3=1-(squeeze(sum(cal3(3,:,2:end))));
x4=1-(squeeze(sum(cal4(3,:,2:end))));

YL=2;
YS=0;

r=0.75;

h1=figure(1);
set(gcf,'Resize','off')
set(h1,'position',[100 100 1000 600]);
dx=0.43;
dy=0.38;
ddx=(1.00-2*dx)/3;
ddy=(0.92-2*dy)/2;

s1=subplot('position',[0*dx+1*ddx+0.01 1*dy+2.5*ddy dx dy]);
plot(x1',c1')
hold on
line(x1,repmat(1,1,9),'Color','black')
title('RCP2.6')
xlim([0 1]);
ylim([YS YL]);
ylabel('Interval ratio')
xlabel('Frequency')
hold off

s2=subplot('position',[1*dx+2*ddx+0.02 1*dy+2.5*ddy dx dy]);
plot(x2',c2')
hold on
line(x2,repmat(1,1,9),'Color','black')
title('RCP4.5')
xlim([0 1]);
ylim([YS YL]);
ylabel('Interval ratio')
xlabel('Frequency')
hold off

s3=subplot('position',[0*dx+1*ddx+0.01 0*dy+1.*ddy dx dy]);
plot(x3',c3')
hold on
line(x3,repmat(1,1,9),'Color','black')
title('RCP6.0')
xlim([0 1]);
ylim([YS YL]);
ylabel('Interval ratio')
xlabel('Frequency')
hold off

s4=subplot('position',[1*dx+2*ddx+0.02 0*dy+1.*ddy dx dy]);
plot(x4',c4')
hold on
line(x4,repmat(1,1,9),'Color','black')
title('RCP8.5')
xlim([0 1]);
ylim([YS YL]);
ylabel('Interval ratio')
xlabel('Frequency')
%legend({'\gamma_u/\delta_G','\gamma_d/\delta_G'})

hold off



fname=['PLOTS/cal_ratio_runmean_20'];
%print(h1,'-dpng','-r300',fname);




%-----------------------------------------------------------------------------
% calculate mean and STD time series
%-----------------------------------------------------------------------------

byear=2020;
eyear=2099;

nmean=eyear-byear+1; %eyear-byear+1;

w1=repmat(ncread(filo1,'weight'),[1 nmean])';
w2=repmat(ncread(filo2,'weight'),[1 nmean])';
w3=repmat(ncread(filo3,'weight'),[1 nmean])';
w4=repmat(ncread(filo4,'weight'),[1 nmean])';

byear=byear-1948-syear+1;
eyear=eyear-1948-syear+1;

p1=sum(w1.*f1((syear+byear):(syear+eyear),:),2);
p2=sum(w2.*f2((syear+byear):(syear+eyear),:),2);
p3=sum(w3.*f3((syear+byear):(syear+eyear),:),2);
p4=sum(w4.*f4((syear+byear):(syear+eyear),:),2);

sigma1=(sqrt(sum(w1.*(f1((syear+byear):(syear+eyear),:)-repmat(p1,[1 nmodel1])).^2,2)));
sigma2=(sqrt(sum(w2.*(f2((syear+byear):(syear+eyear),:)-repmat(p2,[1 nmodel2])).^2,2)));
sigma3=(sqrt(sum(w3.*(f3((syear+byear):(syear+eyear),:)-repmat(p3,[1 nmodel3])).^2,2)));
sigma4=(sqrt(sum(w4.*(f4((syear+byear):(syear+eyear),:)-repmat(p4,[1 nmodel4])).^2,2)));

p1_avg=mean(f1((syear+byear):(syear+eyear),:),2);
p2_avg=mean(f2((syear+byear):(syear+eyear),:),2);
p3_avg=mean(f3((syear+byear):(syear+eyear),:),2);
p4_avg=mean(f4((syear+byear):(syear+eyear),:),2);

sigma1_avg=(sqrt(mean((f1((syear+byear):(syear+eyear),:)-repmat(p1_avg,[1 nmodel1])).^2,2)));
sigma2_avg=(sqrt(mean((f2((syear+byear):(syear+eyear),:)-repmat(p2_avg,[1 nmodel2])).^2,2)));
sigma3_avg=(sqrt(mean((f3((syear+byear):(syear+eyear),:)-repmat(p3_avg,[1 nmodel3])).^2,2)));
sigma4_avg=(sqrt(mean((f4((syear+byear):(syear+eyear),:)-repmat(p4_avg,[1 nmodel4])).^2,2)));


%c1=((sigma1)./(sigma1_avg));
%c2=((sigma2)./(sigma2_avg));
%c3=((sigma3)./(sigma3_avg));
%c4=((sigma4)./(sigma4_avg));

nyear=eyear-byear+1;
nmean=1;

c1=(((squeeze(cal1(1,1,2))+squeeze(cal1(1,2,2)))*sigma1)./((2*icdf('norm',1-squeeze(sum(cal1(3,:,2)))/2,0,1))*sigma1_avg));
c2=(((squeeze(cal2(1,1,2))+squeeze(cal2(1,2,2)))*sigma2)./((2*icdf('norm',1-squeeze(sum(cal2(3,:,2)))/2,0,1))*sigma2_avg));
c3=(((squeeze(cal3(1,1,2))+squeeze(cal3(1,2,2)))*sigma3)./((2*icdf('norm',1-squeeze(sum(cal3(3,:,2)))/2,0,1))*sigma3_avg));
c4=(((squeeze(cal4(1,1,2))+squeeze(cal4(1,2,2)))*sigma4)./((2*icdf('norm',1-squeeze(sum(cal4(3,:,2)))/2,0,1))*sigma4_avg));

C1=mean(reshape(c1,[nmean nyear/nmean]),1)';
C2=mean(reshape(c2,[nmean nyear/nmean]),1)';
C3=mean(reshape(c3,[nmean nyear/nmean]),1)';
C4=mean(reshape(c4,[nmean nyear/nmean]),1)';

byear=2020;
eyear=2099;

h2=figure(2);
plot((byear:nmean:eyear)',[C1 C2 C3 C4],'LineWidth',2)
xlabel('Year')
ylabel('0.9 significance level range ratio')
set(gca,'Ylim',[0 2])
set(gca,'Xlim',[byear-1 eyear+1])
set(gca,'FontSize',12)
hold on
line((byear:nmean:eyear)',repmat(1,1,length((byear:nmean:eyear)')),'Color','black','LineWidth',2)
grid on
l=legend('RCP2.6','RCP4.5','RCP6.0','RCP8.5');
set(l,'FontSize',12)
set(l,'Location','NorthWest')
%set(gca,'Ytick',YS:YD:YL)
hold off




fname=['PLOTS/cal_std_ratio_runmean_20'];
print(h2,'-dpng','-r300',fname);


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

  switch mtype
    case 1
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
  
    case 2
     p1=mean(f1((syear+byear):(syear+eyear),:),1);
     p2=mean(f2((syear+byear):(syear+eyear),:),1);
     p3=mean(f3((syear+byear):(syear+eyear),:),1);
     p4=mean(f4((syear+byear):(syear+eyear),:),1);
  
     sigma1=sqrt(sum(w1(1,:).*(p1-repmat(mean(p1),[1 length(p1)])).^2));
     sigma2=sqrt(sum(w2(1,:).*(p2-repmat(mean(p2),[1 length(p2)])).^2));
     sigma3=sqrt(sum(w3(1,:).*(p3-repmat(mean(p3),[1 length(p3)])).^2));
     sigma4=sqrt(sum(w4(1,:).*(p4-repmat(mean(p4),[1 length(p4)])).^2));
  
     p1_avg=mean(f1((syear+byear):(syear+eyear),:),1);
     p2_avg=mean(f2((syear+byear):(syear+eyear),:),1);
     p3_avg=mean(f3((syear+byear):(syear+eyear),:),1);
     p4_avg=mean(f4((syear+byear):(syear+eyear),:),1);
  
     sigma1_avg=sqrt(mean((p1_avg-repmat(mean(p1_avg),[1 length(p1_avg)])).^2));
     sigma2_avg=sqrt(mean((p2_avg-repmat(mean(p2_avg),[1 length(p2_avg)])).^2));
     sigma3_avg=sqrt(mean((p3_avg-repmat(mean(p3_avg),[1 length(p3_avg)])).^2));
     sigma4_avg=sqrt(mean((p4_avg-repmat(mean(p4_avg),[1 length(p4_avg)])).^2));

    case 3
     p1=sum(w1.*f1((syear+byear):(syear+eyear),:),2);
     p2=sum(w2.*f2((syear+byear):(syear+eyear),:),2);
     p3=sum(w3.*f3((syear+byear):(syear+eyear),:),2);
     p4=sum(w4.*f4((syear+byear):(syear+eyear),:),2);

     sigma1=sqrt(mean(sum(w1.*(f1((syear+byear):(syear+eyear),:)-repmat(p1,[1 nmodel1])).^2,2)));
     sigma2=sqrt(mean(sum(w2.*(f2((syear+byear):(syear+eyear),:)-repmat(p2,[1 nmodel2])).^2,2)));
     sigma3=sqrt(mean(sum(w3.*(f3((syear+byear):(syear+eyear),:)-repmat(p3,[1 nmodel3])).^2,2)));
     sigma4=sqrt(mean(sum(w4.*(f4((syear+byear):(syear+eyear),:)-repmat(p4,[1 nmodel4])).^2,2)));

     p1_avg=mean(f1((syear+byear):(syear+eyear),:),2);
     p2_avg=mean(f2((syear+byear):(syear+eyear),:),2);
     p3_avg=mean(f3((syear+byear):(syear+eyear),:),2);
     p4_avg=mean(f4((syear+byear):(syear+eyear),:),2);

     sigma1_avg=sqrt(mean(mean((f1((syear+byear):(syear+eyear),:)-repmat(p1_avg,[1 nmodel1])).^2,2)));
     sigma2_avg=sqrt(mean(mean((f2((syear+byear):(syear+eyear),:)-repmat(p2_avg,[1 nmodel2])).^2,2)));
     sigma3_avg=sqrt(mean(mean((f3((syear+byear):(syear+eyear),:)-repmat(p3_avg,[1 nmodel3])).^2,2)));
     sigma4_avg=sqrt(mean(mean((f4((syear+byear):(syear+eyear),:)-repmat(p4_avg,[1 nmodel4])).^2,2)));  
  end

c1=(((squeeze(cal1(1,1,2:end))+squeeze(cal1(1,2,2:end)))*sigma1)./((2*icdf('norm',1-squeeze(sum(cal1(3,:,2:end)))/2,0,1))*sigma1_avg))
c2=(((squeeze(cal2(1,1,2:end))+squeeze(cal2(1,2,2:end)))*sigma2)./((2*icdf('norm',1-squeeze(sum(cal2(3,:,2:end)))/2,0,1))*sigma2_avg))
c3=(((squeeze(cal3(1,1,2:end))+squeeze(cal3(1,2,2:end)))*sigma3)./((2*icdf('norm',1-squeeze(sum(cal3(3,:,2:end)))/2,0,1))*sigma3_avg))
c4=(((squeeze(cal4(1,1,2:end))+squeeze(cal4(1,2,2:end)))*sigma4)./((2*icdf('norm',1-squeeze(sum(cal4(3,:,2:end)))/2,0,1))*sigma4_avg))

S1=repmat(sigma1./sigma1_avg,9,1)
S2=repmat(sigma2./sigma2_avg,9,1)
S3=repmat(sigma3./sigma3_avg,9,1)
S4=repmat(sigma4./sigma4_avg,9,1)

g1=(((squeeze(cal1(1,1,2:end))+squeeze(cal1(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal1(3,:,2:end)))/2,0,1))))
g2=(((squeeze(cal2(1,1,2:end))+squeeze(cal2(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal2(3,:,2:end)))/2,0,1))))
g3=(((squeeze(cal3(1,1,2:end))+squeeze(cal3(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal3(3,:,2:end)))/2,0,1))))
g4=(((squeeze(cal4(1,1,2:end))+squeeze(cal4(1,2,2:end))))./((2*icdf('norm',1-squeeze(sum(cal4(3,:,2:end)))/2,0,1))))


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

s1=subplot('position',[0*dx+1*ddx+0.01 1*dy+2.5*ddy dx dy]);
plot(x1',[g1' ; S1' ; c1'],'LineWidth',2)
set(gca,'FontSize',12)
hold on
line(x1,repmat(1,1,9),'Color','black','LineWidth',2)
title('RCP2.6')
xlim([0 1]);
ylim([YS YL]);
ylabel('Ratio')
xlabel('Confidence level')
l=legend('(\gamma_u+\gamma_d)/(2\delta_G)','STD ratio','Interval ratio');
set(l,'FontSize',12)
set(l,'Location','NorthWest')
grid on
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

s2=subplot('position',[1*dx+2*ddx+0.02 1*dy+2.5*ddy dx dy]);
plot(x2',[g2' ; S2' ; c2'],'LineWidth',2)
set(gca,'FontSize',12)
hold on
line(x2,repmat(1,1,9),'Color','black','LineWidth',2)
title('RCP4.5')
xlim([0 1]);
ylim([YS YL]);
ylabel('Ratio')
xlabel('Confidence level')
grid on
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

s3=subplot('position',[0*dx+1*ddx+0.01 0*dy+1.*ddy dx dy]);
plot(x3',[g3' ; S3' ; c3'],'LineWidth',2)
set(gca,'FontSize',12)
hold on
line(x3,repmat(1,1,9),'Color','black','LineWidth',2)
title('RCP6.0')
xlim([0 1]);
ylim([YS YL]);
ylabel('Ratio')
xlabel('Confidence level')
grid on
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

s4=subplot('position',[1*dx+2*ddx+0.02 0*dy+1.*ddy dx dy]);
plot(x4',[g4' ; S4' ; c4'],'LineWidth',2)
set(gca,'FontSize',12)
hold on
line(x4,repmat(1,1,9),'Color','black','LineWidth',2)
title('RCP8.5')
xlim([0 1]);
ylim([YS YL]);
ylabel('Ratio')
xlabel('Confidence level')
grid on
set(gca,'Ytick',YS:YD:YL)
set(gca,'Xtick',0:0.1:1)
hold off

fname=['PLOTS/cal_ratio_all_runmean_20'];
print(h3,'-dpng','-r300',fname);

