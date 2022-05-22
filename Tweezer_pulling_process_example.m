clear;
raw=importdata('B01_tx_mT.dat',' ',5);%import data from traces
calb=importdata('B01_Stiffness.dat',' ',5);%%import calibration data
%%QPD position signal. Note that the x, y of stage are swapped compared to the
%%function generator and PSD calibration
%%Use the stage axis of the calibration and function generator as reference
r=0.585/2;%%radius of bead, unit:um
x_qpd=raw.data (:,1);
y_qpd=raw.data(:,2);
z_qpd=raw.data(:,3);
%% import stage position (unit voltage)   
x_stage=(raw.data(:,5))*3;%%unit micron,middle point as reference (~anchorpoint)
y_stage=(raw.data(:,4))*3;
z_stage=(raw.data(:,6))*1;
%%stiffness of three dimensions (unit: pN/nm)
kappa_x=calb.data(1,2);
kappa_y=calb.data(1,3);
kappa_z=calb.data(1,4);
%%Sensitivity factor of QPD (nm/V)
beta_x=calb.data(1,9)*1000;
beta_y=calb.data(1,10)*1000;
beta_z=calb.data(1,11)*1000;

zi=0.30+r;%%initial trap center distance from the surface ~250-350 nm
time=0.001*[1:length(x_qpd)];%%time unit: sec

t1=146.4*1000;
t2=212*1000;%%define the range of interest by time point
stage_x_trace=x_stage(t1:t2)-17.42;
x_qpd_trace=x_qpd(t1:t2);
stage_y_trace=y_stage(t1:t2)-16.45;
y_qpd_trace=y_qpd(t1:t2);
z_qpd_trace=z_qpd(t1:t2);

%%curve fit to find the offset along the oscillation axis and the qpd
%%voltage offset
tfi=1; %fitting range initial time
tfe=151.37*1000-t1; %fitting range end time

yqfit=y_qpd_trace(tfi:tfe);
y_stage_fit=stage_y_trace(tfi:tfe);

[xData, yData] = prepareCurveData( y_stage_fit, yqfit );

% Set up fittype and options.
ft = fittype( 'b0+b1*(x-a0)+b3*(x-a0).^3+b5*(x-a0).^5+b7*(x-a0).^7', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



%%


xqpd0=0.1785;%%offset of xqpd
yqpd0=fitresult.b0;%%offset of yqpd_0.2677 based on average

x0=0;%%offset of anchorpoint at the x direction unit:um

force_x=kappa_x*beta_x*(x_qpd_trace-xqpd0);%%force on x direction
force_y=kappa_y*beta_y*(y_qpd_trace-yqpd0);

x_bd=(stage_x_trace-x0)*1000+force_x/kappa_x;%%bead position of x in nm
y0=fitresult.a0;%%trap center to anchorpoint y distance based on parity, unit: um

y_bd=(stage_y_trace-y0)*1000+force_y/kappa_y;%% bead y position based on force balance geometry

zv0=8.6922;%%z-qpd voltage where z force=0; i.e. deltaz=(zv0-zqpd)beta_z. calcuated based on fitting the geometry
ztrap=300+r*1000;%% trap to surface distance estimated from calibration, based on fitting geometry:=137.3 with large uncertainty, unit: nm
z_bd=ztrap-beta_z*(zv0-z_qpd_trace);%%bead center z height
force_z=kappa_z*beta_z*(zv0-z_qpd_trace);

lex=(x_bd.^2+y_bd.^2+z_bd.^2).^(1/2)-r*1000;
ftot=((force_x.^2+force_y.^2+force_z.^2).^(1/2));

smoothl=medfilt1(lex,11);
smoothf=medfilt1(ftot,11);

fdf=smoothf(209*1000-t1:211.25*1000-t1);
ldf=smoothl(209*1000-t1:211.25*1000-t1);
%fdf=smoothf(190000-t1:191.54*1000-t1);
%ldf=smoothl(190*1000-t1:191.54*1000-t1);

pfit0=[10 3000 1000];%%initial guess of persistent lenth lp nm, contour length lc nm,  elastic moduls K0 pN
pfitlb=[5 2500 400];%%lower bound of fit
pfitub=[120 3200 5000];%%upper bound of fit
f=@(p) 4.1./p(1)*(0.25./(1-ldf./p(2)+fdf./p(3)).^2-0.25+ldf./p(2)-fdf./p(3));%%trial force based on input parameters guess
obj_F=@(p) (4.1./p(1)*(0.25./(1-ldf./p(2)+fdf./p(3)).^2-0.25+ldf./p(2)-fdf./p(3))-fdf);%%objective function for fitting: difference between experimental force and force based on input guess
%options=optimset('TolFun',1e-10,'TolX',1e-10);
p=lsqnonlin(obj_F,pfit0,pfitlb,pfitub);

ffit=f(p);%%force output from WLC fit

%%calculate the DNA length under force at each point based on WLC model
lp_dna=p(1);%%persistence length of DNA based on fitting
lc_dna=p(2);%%contour length of DNA based on fitting
K0=p(3);%%stretching modulus
r_d=zeros(size(smoothf));
for i=1:size(smoothf,1)
    coeff=[4 4*lp_dna/4.1*smoothf(i)-3 0 -1];%%Coefficient of WLC model, converting to a cubic equation of 1-x/lc_dna+f/K0
    
    r_all=roots(coeff);
    r_d(i)=r_all(r_all>=0);
end
ldna=(1-r_d+smoothf/K0)*lc_dna;%%convert back to extension in nm
lnet=medfilt1(lex-ldna,11);
%plot(lex(120:end),ftot(120:end))

plot(smoothl,smoothf)
hold on;

f_wlc=[0.1:0.01:40]';%%force for plotting WLC fit
rwlc_d=zeros(size(f_wlc));

for j=1:size(f_wlc,1)
    coeff=[4 4*lp_dna/4.11*f_wlc(j)-3 0 -1];%%Coefficient of WLC model, converting to a cubic equation of 1-x/lc_dna+f/K0
    
    rwlc_all=roots(coeff);
    rwlc_d(j)=rwlc_all(rwlc_all>=0);
end
l_wlc=(1-rwlc_d+f_wlc/K0)*lc_dna;

plot(l_wlc,f_wlc)

