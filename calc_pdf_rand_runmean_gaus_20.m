tic

for i=1:3

years=[2005 2065 2099];
byear=years(i); eyear=byear;

Title=[num2str(byear) '-' num2str(eyear)];

R='NCEP';
SLA='AVG';

syear=70;

byear=byear-1948-syear+1;
eyear=eyear-1948-syear+1;

diro='./CAL/';

filo1=[ diro 'cal_RCP26_' R '_' SLA '_runmean_gaus_20.nc'];
filo2=[ diro 'cal_RCP45_' R '_' SLA '_runmean_gaus_20.nc'];
filo3=[ diro 'cal_RCP60_' R '_' SLA '_runmean_gaus_20.nc'];
filo4=[ diro 'cal_RCP85_' R '_' SLA '_runmean_gaus_20.nc'];

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

nmean=eyear-byear+1;

w1=repmat(ncread(filo1,'weight'),[1 nmean])';
w2=repmat(ncread(filo2,'weight'),[1 nmean])';
w3=repmat(ncread(filo3,'weight'),[1 nmean])';
w4=repmat(ncread(filo4,'weight'),[1 nmean])';

%-----------------------------------------------------------------------------
% calculate mean and STD
%-----------------------------------------------------------------------------


p1=sum(w1.*f1((syear+byear):(syear+eyear),:),2);
p2=sum(w2.*f2((syear+byear):(syear+eyear),:),2);
p3=sum(w3.*f3((syear+byear):(syear+eyear),:),2);
p4=sum(w4.*f4((syear+byear):(syear+eyear),:),2);

p1r=sum(w1.*f1(58,:),2);
p2r=sum(w2.*f2(58,:),2);
p3r=sum(w3.*f3(58,:),2);
p4r=sum(w4.*f4(58,:),2);

sigma1=sqrt(sum(w1.*(f1((syear+byear):(syear+eyear),:)-repmat(p1,[1 nmodel1])).^2,2));
sigma2=sqrt(sum(w2.*(f2((syear+byear):(syear+eyear),:)-repmat(p2,[1 nmodel2])).^2,2));
sigma3=sqrt(sum(w3.*(f3((syear+byear):(syear+eyear),:)-repmat(p3,[1 nmodel3])).^2,2));
sigma4=sqrt(sum(w4.*(f4((syear+byear):(syear+eyear),:)-repmat(p4,[1 nmodel4])).^2,2));

x1=zeros(21,nmean);
x2=zeros(21,nmean);
x3=zeros(21,nmean);
x4=zeros(21,nmean);

x1(12:21,:)=repmat(p1',[10 1])+...
           (repmat(flipud(squeeze(cal1(1,1,:))),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(flipud(squeeze(cal1(2,1,:)))*pi/180),[1 nmean])...
           ).*repmat(sigma1',[10 1]);
x1(11,:)=p1;
x1(1:10,:)=repmat(p1',[10 1])-...
           (repmat(squeeze(cal1(1,2,:)),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(squeeze(cal1(2,2,:))*pi/180),[1 nmean])...
           ).*repmat(sigma1',[10 1]);
x1=x1-p1r;

x2(12:21,:)=repmat(p2',[10 1])+...
           (repmat(flipud(squeeze(cal2(1,1,:))),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(flipud(squeeze(cal2(2,1,:)))*pi/180),[1 nmean])...
           ).*repmat(sigma2',[10 1]);
x2(11,:)=p2;
x2(1:10,:)=repmat(p2',[10 1])-...
           (repmat(squeeze(cal2(1,2,:)),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(squeeze(cal2(2,2,:))*pi/180),[1 nmean])...
           ).*repmat(sigma2',[10 1]);
x2=x2-p2r;

x3(12:21,:)=repmat(p3',[10 1])+...
           (repmat(flipud(squeeze(cal3(1,1,:))),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(flipud(squeeze(cal3(2,1,:)))*pi/180),[1 nmean])...
           ).*repmat(sigma3',[10 1]);
x3(11,:)=p3;
x3(1:10,:)=repmat(p3',[10 1])-...
           (repmat(squeeze(cal3(1,2,:)),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(squeeze(cal3(2,2,:))*pi/180),[1 nmean])...
           ).*repmat(sigma3',[10 1]);
x3=x3-p3r;

x4(12:21,:)=repmat(p4',[10 1])+...
           (repmat(flipud(squeeze(cal4(1,1,:))),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(flipud(squeeze(cal4(2,1,:)))*pi/180),[1 nmean])...
           ).*repmat(sigma4',[10 1]);
x4(11,:)=p4;
x4(1:10,:)=repmat(p4',[10 1])-...
           (repmat(squeeze(cal4(1,2,:)),[1 nmean])+...
            repmat(((syear+byear):(syear+eyear)),[10 1]).*...
            repmat(tan(squeeze(cal4(2,2,:))*pi/180),[1 nmean])...
           ).*repmat(sigma4',[10 1]);
x4=x4-p4r;

y1(12:21)=1-flipud(squeeze(cal1(3,1,:)));
y1(11)=0.5;
y1(1:10)=squeeze(cal1(3,2,:));

y2(12:21)=1-flipud(squeeze(cal2(3,1,:)));
y2(11)=0.5;
y2(1:10)=squeeze(cal2(3,2,:));

y3(12:21)=1-flipud(squeeze(cal3(3,1,:)));
y3(11)=0.5;
y3(1:10)=squeeze(cal3(3,2,:));

y4(12:21)=1-flipud(squeeze(cal4(3,1,:)));
y4(11)=0.5;
y4(1:10)=squeeze(cal4(3,2,:));


XI1=x1;
XI2=x2;
XI3=x3;
XI4=x4;

PI1=diff(y1);
PI2=diff(y2);
PI3=diff(y3);
PI4=diff(y4);

n=1e7;
ny=1;
RAND1=zeros(n,ny);
RAND2=zeros(n,ny);
RAND3=zeros(n,ny);
RAND4=zeros(n,ny);

for i=1:ny
   RAND1(:,i)=gen_rand(n,XI1(:,i)',PI1);
   RAND2(:,i)=gen_rand(n,XI2(:,i)',PI2);
   RAND3(:,i)=gen_rand(n,XI3(:,i)',PI3);
   RAND4(:,i)=gen_rand(n,XI4(:,i)',PI4);
end

edges=-0.5:0.05:7.;

[PM1,edges1] = histcounts(mean(RAND1,2),edges, 'Normalization', 'probability');
[PM2,edges2] = histcounts(mean(RAND2,2),edges, 'Normalization', 'probability');
[PM3,edges3] = histcounts(mean(RAND3,2),edges, 'Normalization', 'probability');
[PM4,edges4] = histcounts(mean(RAND4,2),edges, 'Normalization', 'probability');

XM=edges;
DX=diff(edges);
MD=min(XM);
MU=max(XM);

PM1_005=pint(XM,PM1,0.05);
PM2_005=pint(XM,PM2,0.05);
PM3_005=pint(XM,PM3,0.05);
PM4_005=pint(XM,PM4,0.05);

PM1_017=pint(XM,PM1,0.1666);
PM2_017=pint(XM,PM2,0.1666);
PM3_017=pint(XM,PM3,0.1666);
PM4_017=pint(XM,PM4,0.1666);

PM1_05=pint(XM,PM1,0.5);
PM2_05=pint(XM,PM2,0.5);
PM3_05=pint(XM,PM3,0.5);
PM4_05=pint(XM,PM4,0.5);

PM1_083=pint(XM,PM1,0.8333);
PM2_083=pint(XM,PM2,0.8333);
PM3_083=pint(XM,PM3,0.8333);
PM4_083=pint(XM,PM4,0.8333);

PM1_095=pint(XM,PM1,0.95);
PM2_095=pint(XM,PM2,0.95);
PM3_095=pint(XM,PM3,0.95);
PM4_095=pint(XM,PM4,0.95);

Q=[...
[PM1_005 PM1_05 PM1_095] ; ...
[PM2_005 PM2_05 PM2_095] ; ...
[PM3_005 PM3_05 PM3_095] ; ...
[PM4_005 PM4_05 PM4_095]]

toc

s    = [PM1_005 PM2_005 PM3_005 PM4_005 ; ...
         PM1_017 PM2_017 PM3_017 PM4_017 ; ...
         PM1_083 PM2_083 PM3_083 PM4_083 ; ...
         PM1_095 PM2_095 PM3_095 PM4_095];

S2 = [mean([x1(2,:)' p1-clm x1(20,:)'],1) ; ...
      mean([x2(2,:)' p2-clm x2(20,:)'],1) ; ...
      mean([x3(2,:)' p3-clm x3(20,:)'],1) ; ...
      mean([x4(2,:)' p4-clm x4(20,:)'],1)]


%h2=figure(1);
%data = rand(100,4);
%YL=0.3;

fname=['CAL/calc_pdf_rand_RAND_' num2str(byear+2017) '-' num2str(eyear+2017) '_runmean_gaus_20'];
save(fname,'RAND1','RAND2','RAND3','RAND4')

clear RAND1 RAND2 RAND3 RAND4
fname=['CAL/calc_pdf_rand_ALL_' num2str(byear+2017) '-' num2str(eyear+2017) '_runmean_gaus_20'];
save(fname)

clear
end
