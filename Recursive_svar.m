% Gabriel El Mahi, Mariela Huillier, Damien Vliegen 

% Advanced Macroeconomics - Assignment 3
% Recursive SVAR approach

%% Loading data

% To close all the existing figures and clean the workspace and clear the
% terminal command
close all;
clear;
clc;
rng(111); % set seed

% Importing the DATA from our excel worksheet
DATA = readtable('Data.xlsx','Sheet','DATA','Range','A1:H183');

% Remove first entries when we don't have information for tightening variable 
DATA([1:61],:) = [];

% Name variables 
DATA.DATE = table2array(DATA(:,1));
DATA.SIR = table2array(DATA(:,3));
DATA.SPREAD = table2array(DATA(:,4));
DATA.CPI = table2array(DATA(:,5));
DATA.GDP = table2array(DATA(:,6));
DATA.VOL = table2array(DATA(:,7));
DATA.TIGHT = table2array(DATA(:,8));

% Construct variables 
data = table(DATA.DATE(5:end), VariableNames = "DATE");
data.CPI_Y_Y = (DATA.CPI(5:end) - DATA.CPI(1:end-4))./DATA.CPI(1:end-4) * 100;
data.GDP_Y_Y = (DATA.GDP(5:end) - DATA.GDP(1:end-4))./DATA.GDP(1:end-4) * 100;
data.SIR = DATA.SIR(5:end);
data.SPREAD = DATA.SPREAD(5:end);
data.VOL = (DATA.VOL(5:end) - DATA.VOL(1:end-4)) ./ DATA.VOL(1:end-4) * 100;
data.TIGHT = DATA.TIGHT(5:end);
data = table2timetable(data, RowTimes = "DATE");

%% Test number of lags 
[lag, AIC, BIC, HQC]=var_lag([data.GDP_Y_Y data.CPI_Y_Y data.VOL data.TIGHT data.SPREAD data.SIR],10)

%% Plot tightening  
figure()
plot(data.DATE, data.TIGHT)
legend('Lending standards/Tightening')
title({'Evolution of lending standards'})
ylabel('Net percent')
hold off

%% C Matrix identification via Cholesky decomposition, and IRFs

% Selecting the lag
lag = 4;

% Setting the variables (credit shock in 4th position)
Variables = [data.GDP_Y_Y data.CPI_Y_Y data.VOL data.TIGHT data.SPREAD data.SIR];

% OLS Estimation
[Beta,Sigma,ERR,x,T,n] = var_est(Variables,lag);

% C Matrix Identification
C = chol(Sigma, "LOWER");

% Plot IRFs from an increase of loan standards leading to a 0.1 unit increase
% of the spread 
n = 6;
I = eye(n*lag);                         
A = [Beta(2:end,:)'; I(1:end-n,:)]; 

for j = 0:T-1
     P = A^j;
     Phi(:,:,j+1) = P(1:n,1:n);
end
 
for j = 1:size(Phi,3)
     IRF(:,:,j) = Phi(:,:,j)*C;
end

% Create confidence intervals and median
[IRF_lb,IRF_up]=conf_int(Beta,Sigma,T,0.10);
[IRF_med,~]=conf_int(Beta,Sigma,T,0.50);

% Normalize the values (credit shock with impact of +0.1% on credit spread
IRF_lb = IRF_lb/IRF_med(5,4)*0.1;
IRF_up = IRF_up/IRF_med(5,4)*0.1;
IRF_med = IRF_med/IRF_med(5,4)*0.1;

% Plot IRFs and confidence intervals
figure()

subplot(2,3,1)
plotx2(1:T,[squeeze(IRF_lb(1,4,:)) squeeze(IRF_up(1,4,:))]), hold on
plot(squeeze(IRF_med(1,4,:)),'r'), hold on
zero=get(gca,'xlim');
plot(zero,[0 0],'k')
title({'Tightening shock to GDP GROWTH'})
xlabel('Lag')
ylabel('Percent')
axis([0 40 -inf inf])

subplot(2,3,2)
plotx2(1:T,[squeeze(IRF_lb(2,4,:)) squeeze(IRF_up(2,4,:))]), hold on
plot(squeeze(IRF_med(2,4,:)),'r'), hold on
zero=get(gca,'xlim');
plot(zero,[0 0],'k')
title({'Tightening shock to INFLATION'})
xlabel('Lag')
ylabel('Percent')
axis([0 40 -inf inf])

subplot(2,3,3)
plotx2(1:T,[squeeze(IRF_lb(3,4,:)) squeeze(IRF_up(3,4,:))]), hold on
plot(squeeze(IRF_med(3,4,:)),'r'), hold on
zero=get(gca,'xlim');
plot(zero,[0 0],'k')
title({'Tightening shock to CREDIT VOL GROWTH'})
xlabel('Lag')
ylabel('Percent')
axis([0 40 -inf inf])

subplot(2,3,4)
plotx2(1:T,[squeeze(IRF_lb(4,4,:)) squeeze(IRF_up(4,4,:))]), hold on
plot(squeeze(IRF_med(4,4,:)),'r'), hold on
zero=get(gca,'xlim');
plot(zero,[0 0],'k')
title({'Tightening shock to TIGHTENING'})
xlabel('Lag')
ylabel('Percent')
axis([0 40 -inf inf])

subplot(2,3,5)
plotx2(1:T,[squeeze(IRF_lb(5,4,:)) squeeze(IRF_up(5,4,:))]), hold on
plot(squeeze(IRF_med(5,4,:)),'r'), hold on
zero=get(gca,'xlim');
plot(zero,[0 0],'k')
title({'Tightening shock to SPREAD'})
xlabel('Lag')
ylabel('Percent')
axis([0 40 -inf inf])

subplot(2,3,6)
plotx2(1:T,[squeeze(IRF_lb(6,4,:)) squeeze(IRF_up(6,4,:))]), hold on
plot(squeeze(IRF_med(6,4,:)),'r'), hold on
zero=get(gca,'xlim');
plot(zero,[0 0],'k')
title({'Tightening shock to 3M TBILL'})
xlabel('Lag')
ylabel('Percent')
axis([0 40 -inf inf])

%% Historical decomposition

% Compute structural shocks
epsilonhat = ERR*inv(C)';

% Plot (i) historical GDP growth and (ii) GDP growth with TIGHTENING shock only
% Counterfactual GDP (GDP without tightening shocks)

% Method 1
no_str_tight = [epsilonhat(:,1:3) zeros(size(epsilonhat, 1), 1) epsilonhat(:,5:6)];

tight_cf = zeros(T,n);
tight_cf = Variables(1:4,:);

for i = lag + 1 : T
        X = 1;
    for j = 1 : lag
        X = [X tight_cf(i-j,:)];
    end
    tight_cf(i,:) = X*Beta + no_str_tight(i-lag,:)*C';
end

tight_cf = Variables - tight_cf;

figure()
plot(data.DATE, tight_cf(:,1),'--') 
hold on
plot(data.DATE, Variables(:,1))
plot(data.DATE, zeros(T),'--','color','black');
title('Counterfactual Analysis : GDP growth')
legend("Induced by credit tightening shock", "Historical", 'Location','South')
hold off

figure()
plot(data.DATE, tight_cf(:,2),'--') 
hold on
plot(data.DATE, Variables(:,2))
plot(data.DATE, zeros(T),'--','color','black');
title('Counterfactual Analysis : Inflation')
legend("Induced by credit tightening shock", "Historical")
hold off

%% FEVDs

horizon = 40;

% Total MSE
MSE(:,:,1)=C*C';
for h=2:horizon
    MSE(:,:,h)=MSE(:,:,h-1) + Phi(:,:,h)*(C)*C'*Phi(:,:,h)';
end 

% With GDP shock only
MSE_gdp(:,:,1)=C*[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C';
for h=2:horizon
    MSE_gdp(:,:,h)=MSE_gdp(:,:,h-1) + Phi(:,:,h)*C*[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C'*Phi(:,:,h)';
end 

% With INFLATION shock only
MSE_cpi(:,:,1)=C*[0 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C';
for h=2:horizon
    MSE_cpi(:,:,h)=MSE_cpi(:,:,h-1) + Phi(:,:,h)*C*[0 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C'*Phi(:,:,h)';
end 

% With CREDIT VOL shock only
MSE_vol(:,:,1)=C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C';
for h=2:horizon
    MSE_vol(:,:,h)=MSE_vol(:,:,h-1) + Phi(:,:,h)*C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C'*Phi(:,:,h)';
end 

% With TIGHTENING shock only
MSE_tight(:,:,1)=C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C';
for h=2:horizon
    MSE_tight(:,:,h)=MSE_tight(:,:,h-1) + Phi(:,:,h)*C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]*C'*Phi(:,:,h)';
end 

% With SPREAD shock only
MSE_spread(:,:,1)=C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0]*C';
for h=2:horizon
    MSE_spread(:,:,h)=MSE_spread(:,:,h-1) + Phi(:,:,h)*C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0]*C'*Phi(:,:,h)';
end 

% With TBILL shock only
MSE_tbill(:,:,1)=C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1]*C';
for h=2:horizon
    MSE_tbill(:,:,h)=MSE_tbill(:,:,h-1) + Phi(:,:,h)*C*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1]*C'*Phi(:,:,h)';
end 

FEVD = zeros(horizon,n,n);

% Decomposition of variance of GDP Growth
FEVD(:,1,1) = MSE_gdp(1,1,:)./MSE(1,1,:);
FEVD(:,1,2) = MSE_cpi(1,1,:)./MSE(1,1,:);
FEVD(:,1,3) = MSE_vol(1,1,:)./MSE(1,1,:);
FEVD(:,1,4) = MSE_tight(1,1,:)./MSE(1,1,:);
FEVD(:,1,5) = MSE_spread(1,1,:)./MSE(1,1,:);
FEVD(:,1,6) = MSE_tbill(1,1,:)./MSE(1,1,:);

figure()
sp = stackedplot(1:horizon,FEVD(:,1,:)*100);
set(sp, 'DisplayLabels',["GDP growth " "Inflation " "Credit vol growth " "Tightening " "Spread " "3M Tbill "])
title('Forecast Error Variance Decomposition: GDP growth')
xlabel('Horizon')

% Decomposition of variance of Inflation
FEVD(:,2,1) = MSE_gdp(2,2,:)./MSE(2,2,:);
FEVD(:,2,2) = MSE_cpi(2,2,:)./MSE(2,2,:);
FEVD(:,2,3) = MSE_vol(2,2,:)./MSE(2,2,:);
FEVD(:,2,4) = MSE_tight(2,2,:)./MSE(2,2,:);
FEVD(:,2,5) = MSE_spread(2,2,:)./MSE(2,2,:);
FEVD(:,2,6) = MSE_tbill(2,2,:)./MSE(2,2,:);

figure()
ip = stackedplot(1:horizon,FEVD(:,2,:)*100);
set(ip, 'DisplayLabels',["GDP growth " "Inflation " "Credit vol growth " "Tightening " "Spread " "3M Tbill "])
title('Forecast Error Variance Decomposition: Inflation')
xlabel('Horizon')

%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Beta,Sigma,ERR,x,T,n]=var_est(X,lag)

%The function estimates the VAR via OLS

%INPUT: matrix of data, number of lags

%OUTPUT: beta coefficients, variance-covariance matrix, residuals,
%stacked matrix of data

[T,n]=size(X);

%estimate VAR
Yt=X(lag+1:end,:);
x=X(1:end-lag,:);
for l=2:lag
    x=[X(l:end-lag-1+l,:) x];
end
x = [ones(size(x,1),1) x];
Beta=(x'*x)\(x'*Yt);
ERR=Yt-x*Beta;

dem=ERR-sum(ERR,1)./T;
Sigma= (dem'*dem)./(T-1);

% Sigma=cov(ERR);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lag, AIC, BIC, HQC]=var_lag(data,lags)

%The function selects the optimal lag length of a VAR based on Akaike's
%and Bayesian information criteria

%INPUT: matrix of data, desired max num of lags to be consideredin the
%minimization of AIC and BIC

%OUTPUT: highest num of lags among the ones suggested by the mininization
%of AIC, BIC and HQC; values of AIC, BIC, and HQC at different lags

%form matrices for VAR
data=data'; %transpose matrix of data for conformability
[n,T]=size(data);
t=T-lags;
X=zeros(1+n*lags,t);
X(1,:)=ones(1,t);
for j=1:lags
    X(2+n*(j-1):1+n*j,:)=data(:,lags+1-j:end-j);
end

%vectors of results from information criteria (at different lags)
AIC=zeros(lags,1);
BIC=zeros(lags,1);
HQC=zeros(lags,1);

%calculate AIC and BIC for different lags
for j=lags:-1:1
    Y=data(:,lags+1:end);
    X=X(1:1+n*j,:);
    beta=Y*X'/(X*X');
    res=Y-beta*X;
    sigma=((res)*(res'))*(1/t);
    AIC(j,1)=log(det(sigma))+(2*j*(n^2))/t;
    BIC(j,1)=log(det(sigma))+(log(t)*j*(n^2))/t;
    HQC(j,1)=log(det(sigma))+(2*log(log(t))*j*(n^2))/t;
end

%choose optimal lag
[~,lag_AIC]=min(AIC);
[~,lag_BIC]=min(BIC);
[~,lag_HQC]=min(HQC);
disp('Lag selection criteria');
disp('     AIC   BIC   HQC');
lag_opt=[lag_AIC lag_BIC lag_HQC];
disp(lag_opt);

%for convenience, the minimum of the three is automatically chosen
lag=min(lag_opt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PHI_lb,PHI_up]=conf_int(Beta,Sigma,T,alpha)

%The function simulates multiple VARs and computes the relative IRFs to
%create confidence intervals

%INPUT: VAR coefficients, variance-covariance matrix, length of the data,
%significance level

%OUTPUT: lower-bound and upper-bound intervals for the IRFs

%initialise
N = size(Beta,2);
mu = Beta(1,:);
ARlags = size(Beta(2:end,:),1)/N;
R = chol(Sigma,'lower');

for k=1:1000
    
%create sample from initial conditions
X=NaN(200+T,size(Beta,1));
X(1,:)=[1 repmat(mu,1,ARlags)];

err=randn(200+T,N)*R';

for i=2:200+T
    Ynew=X(i-1,:)*Beta+err(i,:);
    X(i,:)=[1 Ynew X(i-1,2:end-N)];
end

Y=X(end-T+1:end,2:N+1);

%estimate VAR
[Betanew,Sigmanew,~,~]=var_est(Y,ARlags);

%invert VAR (Wold decomposition)
[PHInew,~]=var_inv(Y,Betanew,ARlags);

J=size(PHInew,3);

%estimate IRF via Cholesky decomposition
Bo=chol(Sigmanew,'lower');
for j=1:J
    Bnew(:,:,j,k) = PHInew(:,:,j)*Bo;
end

end

Bnew = sort(Bnew,4,'ascend');

PHI_lb=NaN(N,N,J);
PHI_up=NaN(N,N,J);

%define the confidence intervals
for i=1:N
    for h=1:N
        for j=1:J
            PHI_lb(i,h,j)= Bnew(i,h,j,round(alpha*500));
            PHI_up(i,h,j)= Bnew(i,h,j,round(1000-alpha*500));
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PHI,A,P]=var_inv(X,Beta,lag)

%The function inverts the VAR and retrieves the Wold decomposition

%INPUT: matrix of data, matrix of coefficients, number of lags

%OUTPUT: reduced-form MA matrices, companion matrix

[T,n]=size(X);

%create companion matrix
I=eye(n*lag);
A=[Beta(2:end,:)'; I(1:end-n,:)];

%invert companion matrix
for j=0:T-1
    P = A^j;
    PHI(:,:,j+1)=P(1:n,1:n);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot IRFs with given presentation
function plotx2(t,y)
set(gcf,'DefaultAxesColorOrder',[0.8 0.1 0.1;1 0 0;1 0 0;0 0 1]);
cu=y(:,1);
cl=y(:,2);

h=t;
h=h';
hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
set(hh,'edgecolor',[0.5273 0.8047 0.9180]);
set(hh,'facecolor',[0.5273 0.8047 0.9180]);

axis tight
% hold on
% plot(h,y(:,1),'LineWidth',2);

% hold on;
% zz=zeros(size(y,1),1);
% plot(h,zz,'b-');
end