%%% analysis script for Budget and Advec Decomp 
%%% following the Overleaf doc as of May 6

%please start with the relevant data loaded in!
folder = 'weddell60';
saveon = true;% change to folder you want saved in
%% 6-year mean budget 
len = 438;
k=438;
cv = 3.7843e-4;

load('LayerwiseFIXfull.mat')
weights = (Area1 + Area2 + Area3)/sum(Area1 + Area2 + Area3);
load('LayerwiseVHfix438.mat')
phi0 = sum(Ntot.*repmat(weights,[438 1]),2);

load('WedVol438may12.mat')
voladv = advTot; 
load('WeddellFULL5dayFIX.mat')
load('WedSeaIceSnapfix1.mat')
load('Wed61Wk1to438ASSUMP.mat')
advSm = zeros(1, floor(len/k));
bioSm = zeros(1, floor(len/k));
tendSm = zeros(1, floor(len/k));
mixSm = zeros(1, floor(len/k));
resSm = zeros(1, floor(len/k));
corrSm = zeros(1, floor(len/k));
surfSm = zeros(1, floor(len/k));
dilutSm = zeros(1, floor(len/k));
VOLSm = zeros(1, floor(len/k));
PERSm = zeros(1, floor(len/k));
SIexSm = zeros(1, floor(len/k));
partetaCSm = zeros(1, floor(len/k));
toptendSm = zeros(1, floor(len/k));
topcontSm = zeros(1, floor(len/k));
partCetaSm = zeros(1, floor(len/k));
bioetaSm = zeros(1, floor(len/k));
advetaSm = zeros(1, floor(len/k));
volresSm = zeros(1, floor(len/k));
phi0Sm = zeros(1, floor(len/k));
SSHSm = zeros(1, floor(len/k));
epsSISm = zeros(1, floor(len/k));
SIhSm = zeros(1, floor(len/k));


for i= 1
    advSm(:,i) = mean(advTot((i-1)*k+1: i*k));
    resSm(:,i) = mean(resTot((i-1)*k+1: i*k));
    bioSm(:,i) = mean(bioTot((i-1)*k+1: i*k));
    mixSm(:,i) = mean(mixTot((i-1)*k+1: i*k));
    corrSm(:,i) = mean(corrTot((i-1)*k+1: i*k));
    surfSm(:,i) = mean(surfTot((i-1)*k+1: i*k));
    tendSm(:,i) = mean(tendTot((i-1)*k+1: i*k));
    dilutSm(:,i) = mean(dilutTot((i-1)*k+1: i*k));
    VOLSm(:,i) = cv*mean(VOLtot((i-1)*k+1: i*k));
    PERSm(:,i) = mean(FWTot2((i-1)*k+1:i*k));
    SIexSm(:,i) = mean(SIexTot((i-1)*k+1:i*k));
    topcontSm(:,i) = cv*mean(topcont((i-1)*k+1: i*k));
    partetaCSm(:,i) = cv*mean(partetaC((i-1)*k+1: i*k));
    toptendSm(:,i) = cv*mean(toptend((i-1)*k+1: i*k));
    partCetaSm(:,i) = cv*mean(partCeta((i-1)*k+1: i*k));
    advetaSm(:,i) = cv*mean(adveta((i-1)*k+1: i*k));
    bioetaSm(:,i) = cv*mean(bioeta((i-1)*k+1: i*k));
    volresSm(:,i) = mean(Volres((i-1)*k+1: i*k));
    phi0Sm(:,i) = mean(phi0((i-1)*k+1: i*k));
    SSHSm(:,i) = mean(SSHTot((i-1)*k+1: i*k));
    epsSISm(:,i) = mean(-SIres((i-1)*k+1: i*k));
    SIhSm(:,i) = mean(SIhTot((i-1)*k+1: i*k));
end

PERSm = PERSm*cv.*phi0Sm;
SIexSm = SIexSm*cv.*phi0Sm;
epsvol = partetaCSm - corrSm + dilutSm - cv*volresSm.*phi0Sm;
storage = tendSm + partetaCSm - SSHSm*cv.*phi0Sm - SIhSm*cv.*phi0Sm;

resEx = storage - epsvol - epsSISm*cv.*phi0Sm - surfSm - bioSm -(advSm-VOLSm) -SIexSm + PERSm - mixSm 
resALL = resEx+epsSISm*cv.*phi0Sm+epsvol;
% I want a bar graph!

X = [storage (advSm-VOLSm) -PERSm  SIexSm bioSm surfSm resALL];
figure(77)
bar(X)
set(gca, 'XTickLabel', {'Storage', 'vol-con adv', 'PER flow', 'SI flow', 'biology' ,'air-sea', 'residual'})
title('6-year mean DIC budget')
if (saveon)
saveas(gcf, strcat('figures/',folder,'/seasonalbudget.png'))
end

%% Volume conserved annual budget

len = 438;
k=73;
cv = 3.7843e-4;

load('LayerwiseFIXfull.mat')
weights = (Area1 + Area2 + Area3)/sum(Area1 + Area2 + Area3);
load('LayerwiseVHfix438.mat')
phi0 = sum(Ntot.*repmat(weights,[438 1]),2);

load('WedVol438may12.mat')
voladv= advTot; 
load('WeddellFULL5dayFIX.mat')
load('WedSeaIceSnapfix1.mat')
load('Wed61Wk1to438ASSUMP.mat')
advSm = zeros(1, floor(len/k));
bioSm = zeros(1, floor(len/k));
tendSm = zeros(1, floor(len/k));
mixSm = zeros(1, floor(len/k));
resSm = zeros(1, floor(len/k));
corrSm = zeros(1, floor(len/k));
surfSm = zeros(1, floor(len/k));
dilutSm = zeros(1, floor(len/k));
VOLSm = zeros(1, floor(len/k));
PERSm = zeros(1, floor(len/k));
SIexSm = zeros(1, floor(len/k));
partetaCSm = zeros(1, floor(len/k));
toptendSm = zeros(1, floor(len/k));
topcontSm = zeros(1, floor(len/k));
partCetaSm = zeros(1, floor(len/k));
bioetaSm = zeros(1, floor(len/k));
advetaSm = zeros(1, floor(len/k));
volresSm = zeros(1, floor(len/k));
phi0Sm = zeros(1, floor(len/k));
SSHSm = zeros(1, floor(len/k));
epsSISm = zeros(1, floor(len/k));
SIhSm = zeros(1, floor(len/k));

for i= 1:floor(len/k)
    advSm(:,i) = mean(advTot((i-1)*k+1: i*k));
    resSm(:,i) = mean(resTot((i-1)*k+1: i*k));
    bioSm(:,i) = mean(bioTot((i-1)*k+1: i*k));
    mixSm(:,i) = mean(mixTot((i-1)*k+1: i*k));
    corrSm(:,i) = mean(corrTot((i-1)*k+1: i*k));
    surfSm(:,i) = mean(surfTot((i-1)*k+1: i*k));
    tendSm(:,i) = mean(tendTot((i-1)*k+1: i*k));
    dilutSm(:,i) = mean(dilutTot((i-1)*k+1: i*k));
    VOLSm(:,i) = cv*mean(VOLtot((i-1)*k+1: i*k));
    PERSm(:,i) = mean(FWTot2((i-1)*k+1:i*k));
    SIexSm(:,i) = mean(SIexTot((i-1)*k+1:i*k));
    topcontSm(:,i) = cv*mean(topcont((i-1)*k+1: i*k));
    partetaCSm(:,i) = cv*mean(partetaC((i-1)*k+1: i*k));
    toptendSm(:,i) = cv*mean(toptend((i-1)*k+1: i*k));
    partCetaSm(:,i) = cv*mean(partCeta((i-1)*k+1: i*k));
    advetaSm(:,i) = cv*mean(adveta((i-1)*k+1: i*k));
    bioetaSm(:,i) = cv*mean(bioeta((i-1)*k+1: i*k));
    volresSm(:,i) = mean(Volres((i-1)*k+1: i*k));
    phi0Sm(:,i) = mean(phi0((i-1)*k+1: i*k));
    SSHSm(:,i) = mean(SSHTot((i-1)*k+1: i*k));
    epsSISm(:,i) = mean(-SIres((i-1)*k+1: i*k));
    SIhSm(:,i) = mean(SIhTot((i-1)*k+1: i*k));
end
PERSm = PERSm*cv.*phi0Sm;
SIexSm = SIexSm*cv.*phi0Sm;
epsvol = partetaCSm - corrSm + dilutSm - cv*volresSm.*phi0Sm;

figure(35)
yea = 1:1:len/k;
hold on

storage = tendSm + partetaCSm - SSHSm*cv.*phi0Sm - SIhSm*cv.*phi0Sm;

plot(storage)
%plot(yea,tendSm)
plot(yea,advSm-VOLSm)
plot(yea, -PERSm)
plot(yea, SIexSm)

resEx = storage - epsvol - epsSISm*cv.*phi0Sm - surfSm - bioSm -(advSm-VOLSm) -SIexSm + PERSm - mixSm 

%plot(yea,partetaCSm)
plot(yea,bioSm)
plot(yea,surfSm)

resALL = resEx+epsSISm*cv.*phi0Sm+epsvol;
plot(yea,resALL)
%plot(yea,resEx)
%plot(yea,epsSISm*cv.*phi0Sm)
%plot(yea,epsvol)

legend('Storage', 'Vol-con advection','P-E+R driven flow','SI ex driven flow','Biological S/S', 'Air-sea flux','total res')%,'100*Resid')


set(gca,'xtick',1:6)
title('Weddell Annual Volume-Conserved Budget')%+string(lat))
ylabel('Tg/yr')

xlabel('Years')
if (saveon)
saveas(gcf, strcat('figures/',folder,'/VolConbudget.png'))
end
%%  Seasonality, Volume Conserved 
% using a moving mean with 6 points (monthly smoothing)
% averaging over all 6 years first... 
% then doing the smoothing 






%% advection decomposition ... seasonally smoothed for now! (will need to rewrite, likely)

horSm = zeros(1, seas);
nvhSm = zeros(1, seas);
volSm = zeros(1, seas);
ovrSm = zeros(1, seas);
advSm = zeros(1, seas);

for i= 1:6
    for k = 1:4
    sea = 4*i+k-4;
    t1 = 73*(i-1)+ 18*(k-1)+1 ;
    if (k==4)
        t2 = 73*(i-1)+ 73;
    end
    if (k < 4)
        t2 = 73*(i-1) + 18*(k);
    end
    advSm(:,sea) = mean(advTot(t1:t2));
    horSm(:,sea) = cv*mean(sum(REStot(t1:t2,:),2));
    nvhSm(:,sea) = cv*mean(sum(NVHtot(t1:t2,:),2));
    ovrSm(:,sea) = cv*mean(sum(OVRtot(t1:t2,:),2));
    volSm(:,sea) = cv*mean(VOLtot(t1:t2));
    end
end

figure(303)
yea = 0.125:0.25:5.875;
plot(yea,advSm)
hold on
plot(yea,nvhSm)
plot(yea,volSm)
plot(yea,(nvhSm-volSm))

legend('Budget Adv','Offline Adv', 'Net-vol Adv', 'Vol-con Adv')% 'budget')
ylabel('Tg C/yr')
grid on
if (saveon)
saveas(gcf, strcat('figures/',folder,'/Voldecomp.png'))
end

%% Horiz/Overturning, Seasonally Smoothed
figure(305)
yea = 0.125:0.25:5.875;
plot(yea,nvhSm-volSm)
hold on
plot(yea,ovrSm - volSm)
plot(yea,horSm)

legend('Vol-con Adv', 'Overturning', 'Horizontal')% 'budget')
ylabel('Tg C/yr')
grid on
if (saveon)
saveas(gcf, strcat('figures/',folder,'/reynoldsSmooth.png'))
end

%% Horiz/Overturning Decomp, unsmoothed
figure(304)
years = (1:438)/73;
plot(years,cv*(sum(NVHtot,2)-VOLtot))
hold on
plot(years,cv*(sum(OVRtot,2)-VOLtot))

plot(years,cv*(sum(REStot,2)))
%plot(years,advTot(1:438) - cv*VOLtot)
%plot(advnotfix - cv*sum(NVHtot,2))
legend('Vol-con Adv', 'Overturning', 'Horizontal')% 'budget')
ylabel('Tg C/yr')
grid on
if (saveon)
saveas(gcf, strcat('figures/',folder,'/reynoldsRough.png'))
end

%% residual decomposition
bg= 1;
en = 438;
ye = (1:438)/73;
figure(306)
plot(ye,cv*sum(REStot(bg:en,1:3),2))
hold on
plot(ye,cv*sum(REStot(bg:en,4:5),2))
plot(ye,cv*sum(REStot(bg:en,6:8),2))
plot(ye,cv*sum(REStot(bg:en,9:11),2))
plot(yea,(horSm),'k--')
legend('SW', 'WW', 'CDW','AABW','total')
grid on
ylabel('Tg DIC/ yr')
xlabel('Years')

grid on
title('Layerwise Horizontal Transport into WG')

if (saveon)
saveas(gcf, strcat('figures/',folder,'/horizlayers.png'))
end
