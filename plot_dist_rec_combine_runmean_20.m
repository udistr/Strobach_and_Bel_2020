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
s1=[PM1_05 PM2_05 PM3_05 PM4_05];
h2=figure;
h    = boxplot(data);
set(gca,'ylim',[0 4])
set(gca,'FontSize',FS)
set(gca,'XTickLabel',{'RCP2.6','RCP4.5','RCP6.0','RCP8.5'})
ylabel('GMST change [^{\circ}C]')
set(h(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
set(h(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
set(h(3,:),{'Ydata'},num2cell(s([end end],:),1)')
set(h(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
set(h(5,:),{'Ydata'},num2cell(s([2 3 3 2 2],:),1)')
set(h(6,:),{'Ydata'},num2cell(s1([1 1],:),1)')
title(Title)

fname=['PLOTS/box_year_' num2str(byear+2017) '-' num2str(eyear+2017) '_runmean_20'];
%print(h2,'-dpng','-r300',fname);
savefig(h2,fname)

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
s1=[PM1_05 PM2_05 PM3_05 PM4_05];
h2=figure;
h    = boxplot(data);
set(gca,'ylim',[0 4])
set(gca,'FontSize',FS)
set(gca,'XTickLabel',{'RCP2.6','RCP4.5','RCP6.0','RCP8.5'})
ylabel('GMST change [^{\circ}C]')
set(h(1,:),{'Ydata'},num2cell(s(end-1:end,:),1)')
set(h(2,:),{'Ydata'},num2cell(s(2:-1:1,:),1)')
set(h(3,:),{'Ydata'},num2cell(s([end end],:),1)')
set(h(4,:),{'Ydata'},num2cell(s([1 1],:),1)')
set(h(5,:),{'Ydata'},num2cell(s([2 3 3 2 2],:),1)')
set(h(6,:),{'Ydata'},num2cell(s1([1 1],:),1)')
title(Title)

fname=['PLOTS/box_year_' num2str(byear+2017) '-' num2str(eyear+2017) '_runmean_20'];
%print(h2,'-dpng','-r300',fname);
savefig(h2,fname)

tex_table = latex(vpa(sym([{'RCP26' ; 'RCP45' ; 'RCP60' ; 'RCP85'} num2cell(round([ Q1 Q2],2))])))  
tex_table = latex(vpa(sym([{'RCP26' ; 'RCP45' ; 'RCP60' ; 'RCP85'} num2cell(round([ Q1(:,3)-Q1(:,1) Q2(:,3)-Q2(:,1)],2))])))

load CAL/calc_pdf_rand_ALL_2005-2005_runmean_20.mat
XMcc=XM;
% PM$ was chosen in all case because it uses most of the models (33)
% the other 3 are subsets of PM4
PM1c=PM4;
PM2c=PM4;
PM3c=PM4;
PM4c=PM4;

dr=0.4;
r1=0.05;
r2=0.55;
YL=0.6;
XL=4.5;
%MD=MD-0.5;
%MU=MU-0.5;

h=figure;
set(gcf,'Resize','off')
set(h,'position',[100 100 1000 600]);
dx=0.43;
dy=0.37;
ddx=(1.00-2*dx)/3;
ddy=(0.90-2*dy)/2; 

FS=12;

s1=subplot('position',[0*dx+1*ddx+0.02 1*dy+2.5*ddy dx dy]);
set(gca,'FontSize',FS)
hold on
for ii=1:length(XMaa)-1
  if (PM1a(ii)~=0)
    h1=rectangle('position',[XMaa(ii)+(r1*DX(ii)) 0 dr*DX(ii) PM1a(ii)]);
    h1.FaceColor='red';
    h1.EdgeColor='none';
  end
end
for ii=1:length(XMbb)-1
  if (PM1b(ii)~=0)
    h2=rectangle('position',[XMbb(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM1b(ii)]);
    h2.FaceColor='blue';
    h2.EdgeColor='none';
  end
end
for ii=1:length(XMbb)-1
  if (PM1c(ii)~=0)
    h3=rectangle('position',[XMcc(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM1c(ii)]);
    h3.FaceColor='black';
    h3.EdgeColor='none';
  end
end
title('RCP2.6')
set(gca,'XLim',[MD MU]);
set(gca,'XTick',MD:1:MU);
ylim([0 YL]);
xlim([-0.5 XL]);
hold off

s2=subplot('position',[1*dx+2*ddx+0.02 1*dy+2.5*ddy dx dy]);
set(gca,'FontSize',FS)
hold on
for ii =1:length(XMaa)-1
  if (PM2a(ii)~=0)
    h1=rectangle('position',[XMaa(ii)+(r1*DX(ii)) 0 dr*DX(ii) PM2a(ii)]);
    h1.FaceColor='red';
    h1.EdgeColor='none';
  end
end
for ii=1:length(XMbb)-1
  if (PM2b(ii)~=0)
    h2=rectangle('position',[XMbb(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM2b(ii)]);
    h2.FaceColor='blue';
    h2.EdgeColor='none';
  end
end
for ii=1:length(XMcc)-1
  if (PM1c(ii)~=0)
    h3=rectangle('position',[XMcc(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM2c(ii)]);
    h3.FaceColor='black';
    h3.EdgeColor='none';
  end
end
title('RCP4.5')
set(gca,'XLim',[MD MU]);
set(gca,'XTick',MD:1:MU);
ylim([0 YL]);
xlim([-0.5 XL]);
hold off

s3=subplot('position',[0*dx+1*ddx+0.02 0*dy+1.3*ddy dx dy]);
set(gca,'FontSize',FS)
hold on
for ii=1:length(XMaa)-1
  if (PM3a(ii)~=0)
    h1=rectangle('position',[XMaa(ii)+(r1*DX(ii)) 0 dr*DX(ii) PM3a(ii)]);
    h1.FaceColor='red';
    h1.EdgeColor='none';
  end
end
for ii=1:length(XMbb)-1
  if (PM3b(ii)~=0)
    h2=rectangle('position',[XMbb(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM3b(ii)]);
    h2.FaceColor='blue';
    h2.EdgeColor='none';
  end
end
for ii=1:length(XMcc)-1
  if (PM3c(ii)~=0)
    h3=rectangle('position',[XMcc(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM3c(ii)]);
    h3.FaceColor='black';
    h3.EdgeColor='none';
  end
end
title('RCP6.0')
set(gca,'XLim',[MD MU]);
set(gca,'XTick',MD:1:MU);
ylim([0 YL]);
xlim([-0.5 XL]);
hold off

s4=subplot('position',[1*dx+2*ddx+0.02 0*dy+1.3*ddy dx dy]);
set(gca,'FontSize',FS)
hold on
for ii=1:length(XMaa)-1
  if (PM4a(ii)~=0)
    h1=rectangle('position',[XMaa(ii)+(r1*DX(ii)) 0 dr*DX(ii) PM4a(ii)]);
    h1.FaceColor='red';
    h1.EdgeColor='none';
  end
end
for ii=1:length(XMbb)-1
  if (PM4b(ii)~=0)
    h2=rectangle('position',[XMbb(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM4b(ii)]);
    h2.FaceColor='blue';
    h2.EdgeColor='none';
  end
end
for ii=1:length(XMcc)-1
  if (PM4c(ii)~=0)
    h3=rectangle('position',[XMcc(ii)+(r2*DX(ii)) 0 dr*DX(ii) PM4c(ii)]);
    h3.FaceColor='black';
    h3.EdgeColor='none';
  end
end
title('RCP8.5')
set(gca,'XLim',[MD MU]);
set(gca,'XTick',MD:1:MU);
ylim([0 YL]);
xlim([-0.5 XL]);
hold off

hly=ylabel('Probability');
set(hly,'FontSize',20)
set(hly,'Position',[-6.4405 0.6600 -1])
hlx=xlabel('GMST change [^{\circ}C]');
set(hlx,'FontSize',20)
set(hlx,'Position',[-0.3 -0.05 -1])

fname=['PLOTS/pdf_year_combine_runmean_20'];
print(h,'-dpng','-r300',fname);
savefig(h,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% skewness and kurtosis - annual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sk1a=skew(X1a,Y1a);
sk2a=skew(X2a,Y2a);
sk3a=skew(X3a,Y3a);
sk4a=skew(X4a,Y4a);

sk1b=skew(X1b,Y1b);
sk2b=skew(X2b,Y2b);
sk3b=skew(X3b,Y3b);
sk4b=skew(X4b,Y4b);

kr1a=kurt(X1a,Y1a)-3;
kr2a=kurt(X2a,Y2a)-3;
kr3a=kurt(X3a,Y3a)-3;
kr4a=kurt(X4a,Y4a)-3;

kr1b=kurt(X1b,Y1b)-3;
kr2b=kurt(X2b,Y2b)-3;
kr3b=kurt(X3b,Y3b)-3;
kr4b=kurt(X4b,Y4b)-3;


A1=[[sk1a kr1a sk1b kr1b];...
   [sk2a kr2a sk2b kr2b];...
   [sk3a kr3a sk3b kr3b];...
   [sk4a kr4a sk4b kr4b]]
tex_table = latex(vpa(sym([{'RCP26' ; 'RCP45' ; 'RCP60' ; 'RCP85'} num2cell(round([A1(:,1:2)],2))])))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% skewness and kurtosis - 20year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xm=XM(1:end-1)+diff(XM)/2;
sk1a=skew(xm,PM1a);
sk2a=skew(xm,PM2a);
sk3a=skew(xm,PM3a);
sk4a=skew(xm,PM4a);

sk1b=skew(xm,PM1a);
sk2b=skew(xm,PM2b);
sk3b=skew(xm,PM3b);
sk4b=skew(xm,PM4b);

kr1a=kurt(xm,PM1a)-3;
kr2a=kurt(xm,PM2a)-3;
kr3a=kurt(xm,PM3a)-3;
kr4a=kurt(xm,PM4a)-3;

kr1b=kurt(xm,PM1a)-3;
kr2b=kurt(xm,PM2b)-3;
kr3b=kurt(xm,PM3b)-3;
kr4b=kurt(xm,PM4b)-3;


A2=[[sk1a kr1a sk1b kr1b];...
   [sk2a kr2a sk2b kr2b];...
   [sk3a kr3a sk3b kr3b];...
   [sk4a kr4a sk4b kr4b]]
tex_table = latex(vpa(sym([{'RCP26' ; 'RCP45' ; 'RCP60' ; 'RCP85'} num2cell(round([A2],2))])))








%exit 
