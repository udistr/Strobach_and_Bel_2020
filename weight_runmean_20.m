R='NCEP';
SLA='EGA';

filo1=[ 'CAL/cal_RCP26_' R '_' SLA '_runmean_20.nc'];
filo2=[ 'CAL/cal_RCP45_' R '_' SLA '_runmean_20.nc'];
filo3=[ 'CAL/cal_RCP60_' R '_' SLA '_runmean_20.nc'];
filo4=[ 'CAL/cal_RCP85_' R '_' SLA '_runmean_20.nc'];

w1=ncread(filo1,'weight');
w2=ncread(filo2,'weight');
w3=ncread(filo3,'weight');
w4=ncread(filo4,'weight');

n1=ncread(filo1,'model');
n2=ncread(filo2,'model');
n3=ncread(filo3,'model');
n4=ncread(filo4,'model');

W1=zeros(1,34);
W2=zeros(1,34);
W3=zeros(1,34);
W4=zeros(1,34);

W1(n1+1)=w1;
W2(n2+1)=w2;
W3(n3+1)=w3;
W4(n4+1)=w4;

models={'ACCESS1-0','ACCESS1-3','CCSM4','CMCC-CESM','CMCC-CMS','CMCC-CM','CNRM-CM5', ...
        'CSIRO-Mk3-6-0','CanESM2','GISS-E2-H-CC','GISS-E2-H','GISS-E2-R-CC', ...
        'GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR', ...
        'IPSL-CM5B-LR','MIROC-ESM-CHEM','MIROC-ESM','MIROC5','MPI-ESM-LR','MPI-ESM-MR', ...
        'MRI-CGCM3','MRI-ESM1','NorESM1-ME','NorESM1-M','CESM1-BGC','CESM1-CAM5', ...
        'EC-EARTH','FGOALS-g2','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M'}


W=vpa(sym(round([W1' W2' W3' W4'],2)));
W(W==0)=str2sym('X')

tex_table = latex([str2sym(models') W])




