% Advanced Macroeconomics I - Assignment 3
% Credit Market Shock and Real Activity
% Method 2: Instrumental Variable Approach

%% Part I: Instrument = "News" (following Mumtaz et al. (2018))

% To close all the existing figures and clean the workspace and clear the
% terminal command
close all;
clear;
clc;

%% Data and Parameters Setting %%

% Importing the DATA from excel worksheet
DATA = readtable('Data.xlsx','Sheet','DATA','Range','A1:M183');

% Definition of vectors of variables from the DATA
DATE = table2array(DATA(:,1));
SIR = table2array(DATA(:,3));
SPREAD = table2array(DATA(:,4));
CPI = table2array(DATA(:,5));
GDP = table2array(DATA(:,6));
VOL = table2array(DATA(:,7));
TIGHT = table2array(DATA(:,8));
NEWS = table2array(DATA(:,9));

figure
plot(DATE,NEWS)
hold on
yyaxis right
plot(DATE,TIGHT)
title("Instruments Considered")
legend("News (Left axis)","Credit Tightening (Right axis - %)")
hold off

% Defining the time interval of interest
Start = 21;
Stop = 152;

% Calculating the AR model variables of interest
Growth = (GDP(Start:Stop) - GDP(Start-4:Stop-4)) ./ GDP(Start-4:Stop-4) * 100;
Inflation = (CPI(Start:Stop) - CPI(Start-4:Stop-4)) ./ CPI(Start-4:Stop-4) * 100;
Rate = SIR(Start:Stop);
Spread = SPREAD(Start:Stop);
Vol = (VOL(Start:Stop) - VOL(Start-4:Stop-4)) ./ VOL(Start-4:Stop-4) * 100;
Variables = [Growth Inflation Rate Spread Vol];
K = size(Variables,2);

% Potential Instruments
Tight = TIGHT(Start:Stop);
News = NEWS(Start:Stop);

% Selecting the instrument
Instrument = News; % Textual measure of credit supply shock

% Select variable instrumented
inst = 4; % Credit spread

% Defining the lag
% [lag, AIC, BIC, HQC] = var_lag(Variables,8);
lag = 4; % selected on the basis of Mumtaz et al. (2018)

% VAR Estimation (OLS)
[BETA,SIGMA,ERR,X] = var_est(Variables,lag);

%% Identification of C matrix %%

% Identifying the relevant column of the "C" matrix (converts credit supply shock 
% into reduced form shocks)
Gamma = (ERR'*Instrument(lag+1:end))/length(Instrument(lag+1:end));

% Normalizing the relevant column of the "C" matrix (shock leads to 1 unit 
% increase of instrumented variable)
C_inst = Gamma / Gamma(inst,1);

%% Simple IRF %%
horizon = 40;

% Building the companion matrix
A = zeros(K*lag);
A(1:K,:) = BETA(2:end,:)';
A = A + diag(ones((lag-1)*K,1),-K);

% Calculation
IRF = zeros(horizon+1,K);
Shock = zeros(K*lag,1);
Shock(1:K,1) = C_inst;
for i = 1 : horizon+1
    response = Shock'*(A^(i-1))';
    IRF(i,:) = response(1,1:K);
end

figure
subplot(2,3,1)
hold on
plot(0:1:horizon,IRF(1:horizon+1,1))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("GDP Growth Response to Credit Shock")
hold off
subplot(2,3,2)
hold on
plot(0:1:horizon,IRF(1:horizon+1,2))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Inflation Response to Credit Shock")
hold off
subplot(2,3,3)
hold on
plot(0:1:horizon,IRF(1:horizon+1,3))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Policy Rate Response to Credit Shock")
hold off
subplot(2,3,4)
hold on
plot(0:1:horizon,IRF(1:horizon+1,4))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Rate Spread Response to Credit Shock")
hold off
subplot(2,3,5)
hold on
plot(0:1:horizon,IRF(1:horizon+1,5))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Volume Growth Response to Credit Shock")
hold off

%% IRF with confidence intervals (bootstrap) %%

% Bootstrap
replications = 1000;
alpha = 0.1;
vector = 1:1:size(ERR,1);
reshuffled_index = RESHUFFLE(vector',10*replications);
ERR_RESHUFFLE = zeros(size(ERR,1),size(ERR,2),replications);
Instrument_RESHUFFLE = zeros(size(Instrument,1),replications);
j = 1;
k = 1;
while j <= 10*replications && k<= replications
    Instrument_RESHUFFLE_TENTATIVE(1:lag,1) = zeros(lag,1);
    for i = 1 : size(ERR,1)
        ERR_RESHUFFLE_TENTATIVE(i,:) = ERR(reshuffled_index(i,1,j),:);
        Instrument_RESHUFFLE_TENTATIVE(i+lag,1) = Instrument(reshuffled_index(i,1,j)+lag); 
    end
    if corr(ERR_RESHUFFLE_TENTATIVE(:,4),Instrument_RESHUFFLE_TENTATIVE(lag+1:end,1)) > 0.3
        ERR_RESHUFFLE(:,:,k) = ERR_RESHUFFLE_TENTATIVE;
        Instrument_RESHUFFLE(:,k) = Instrument_RESHUFFLE_TENTATIVE;
        k = k+1;
    end
    j = j+1;
end

for k = 1 : replications
    Y_star = zeros(size(Variables));
    Y_star(1:lag,:) = Variables(1:lag,:);
    for i = lag+1 : size(Variables,1)
        for j = 1 : lag
            Y_star(i,:) = Y_star(i,:) + Y_star(i-j,:) * BETA(2+(j-1)*K:1+j*K,:);
        end
        Y_star(i,:) = Y_star(i,:) + BETA(1,:) + ERR_RESHUFFLE(i-lag,:,k);
    end
    [BETA_star(:,:,k),SIGMA_star,ERR_star,X_star] = var_est(Y_star,lag);
    C_inst_star(:,k) = (ERR_RESHUFFLE(:,:,k)'*Instrument_RESHUFFLE(lag+1:end,k))/length(Instrument_RESHUFFLE(lag+1:end,k));
    C_inst_star(:,k) = C_inst_star(:,k) / C_inst_star(inst,k);
end

% IRF Calculation
horizon = 40;

IRF_star = zeros(horizon+1,K,replications);
for k = 1 : replications
    % Building the companion matrix
    A = zeros(K*lag);
    A(1:K,:) = BETA_star(2:end,:,k)';
    A = A + diag(ones((lag-1)*K,1),-K);
    % Initialization
    Shock = zeros(K*lag,1);
    Shock(1:K,1) = C_inst_star(:,k);
    % Loop on time steps
    for i = 1 : horizon+1
        response = Shock'*(A^(i-1))';
        IRF_star(i,:,k) = response(1,1:K);
    end
end

impact = 0.1 ; % Impact of initial credit shock on instrumented variable
 
IRF_new = impact * sort(IRF_star,3,'ascend');
IRF_low = IRF_new(:,:,replications * alpha/2);
IRF_med = IRF_new(:,:,replications * 0.5);
IRF_high = IRF_new(:,:,replications * (1-alpha/2));

figure 
subplot(2,3,1)
hold on
plotx2(0:horizon,[IRF_low(:,1) IRF_high(:,1)])
plot(0:1:horizon,IRF_med(:,1),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("GDP Growth Response to Credit Supply Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,2)
hold on
plotx2(0:horizon,[IRF_low(:,2) IRF_high(:,2)])
plot(0:1:horizon,IRF_med(:,2),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Inflation Response to Credit Supply Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,3)
hold on
plotx2(0:horizon,[IRF_low(:,3) IRF_high(:,3)])
plot(0:1:horizon,IRF_med(:,3),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Policy Rate Response to Credit Supply Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,4)
hold on
plotx2(0:horizon,[IRF_low(:,4) IRF_high(:,4)])
plot(0:1:horizon,IRF_med(:,4),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Rate Spread Response to Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,5)
hold on
plotx2(0:horizon,[IRF_low(:,5) IRF_high(:,5)])
plot(0:1:horizon,IRF_med(:,5),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Volume Growth Response to Credit Shock")
legend("Confidence interval (90%)","Median")
hold off

%% Historical Decomposition %%

% Retrieving structural shock on instrumented variable
epsilon_inst = ERR * inv(SIGMA) * Gamma;

% Normalizing the structural shock on instrumented variable
epsilon_inst = epsilon_inst / sqrt(C_inst' * inv(SIGMA) * C_inst);

figure
plot(DATE(Stop-size(epsilon_inst)+1:Stop),epsilon_inst)
title("Structural Credit Shocks (normalized to epsilon/sigma)")

% Counterfactual if only credit supply shocks had occured
Var_CF = zeros(size(Variables,1),size(Variables,2));
for i = lag+1 : size(Var_CF,1)
    X = 0;
    for j = 1 : lag
        X = [X Var_CF(i-j,:)];
    end
    Var_CF(i,:) = X*BETA + ERR(i-lag,:)*inv(SIGMA)*Gamma*Gamma'/(Gamma'*inv(SIGMA_star)*Gamma); 
end

figure
subplot(2,1,1)
plot(DATE(Stop-size(Var_CF,1)+1:Stop),Var_CF(:,1),'--')
hold on
plot(DATE(Stop-size(Variables,1)+1:Stop),Variables(:,1))
plot(DATE(Stop-size(Variables,1)+1:Stop),zeros(length(Variables(:,1))),"--","color","black")
title("GDP Growth (YoY %)")
legend("Credit Supply Shock Induced","Historical")
hold off
subplot(2,1,2)
plot(DATE(Stop-size(Var_CF,1)+1:Stop),Var_CF(:,2),'--')
hold on
plot(DATE(Stop-size(Variables,1)+1:Stop),Variables(:,2))
plot(DATE(Stop-size(Variables,1)+1:Stop),zeros(length(Variables(:,2))),"--","color","black")
title("Inflation (YoY %)")
legend("Credit Supply Shock Induced","Historical")
hold off

%% Forecast Error Variance Decomposition %%

% Parameter setting
horizon = 40;

% Building the companion matrix
I = eye(lag*K);
A = [BETA(2:end,:)';I(1:end-K,:)];

% Phi matrices
Phi = zeros(K,K,horizon+1);
for i = 1 : horizon+1
    matrix = A^(i-1);
    Phi(:,:,i) = matrix(1:K,1:K);
end

% Calculating the MSE
MSE = zeros(K,K,horizon+1);
MSE(:,:,1) = SIGMA;
for i = 1 : horizon
    MSE(:,:,i+1) = MSE(:,:,i+1) + Phi(:,:,i+1) * MSE(:,:,1) * Phi(:,:,i+1)';
end

% Calculating the MSE due to structural credit shock only
MSE_cs = zeros(K,K,horizon+1);
MSE_cs(:,:,1) = Gamma*Gamma'/(Gamma'*inv(SIGMA)*Gamma);
for i = 1 : horizon
    MSE_cs(:,:,i+1) = MSE_cs(:,:,i+1) + Phi(:,:,i+1) * MSE_cs(:,:,1) * Phi(:,:,i+1)';
end

FEVD = zeros(horizon+1,K);
for i = 1 : K
    FEVD(:,i) = MSE_cs(i,i,:) ./ MSE(i,i,:)*100;
end

figure
subplot(2,3,1)
plot(0:horizon,FEVD(:,1))
title("FEVD : Credit Supply Shock on GDP Growth")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,2)
plot(0:horizon,FEVD(:,2))
title("FEVD : Credit Supply Shock on Inflation")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,3)
plot(0:horizon,FEVD(:,3))
title("FEVD : Credit Supply Shock on Policy Rate")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,4)
plot(0:horizon,FEVD(:,4))
title("FEVD : Credit Supply Shock on Credit Spread")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,5)
plot(0:horizon,FEVD(:,5))
title("FEVD : Credit Supply Shock on Credit Volume Growth")
xlabel("Horizon")
ylabel("% MSE")

%% Part II: Instrument = "Tightening" (following Lown and Morgan (2001))

% To close all the existing figures and clean the workspace and clear the
% terminal command
clear;

%% Data and Parameters Setting %%

% Importing the DATA from excel worksheet
DATA = readtable('Data.xlsx','Sheet','DATA','Range','A1:M183');

% Definition of vectors of variables from the DATA
DATE = table2array(DATA(:,1));
SIR = table2array(DATA(:,3));
SPREAD = table2array(DATA(:,4));
CPI = table2array(DATA(:,5));
GDP = table2array(DATA(:,6));
VOL = table2array(DATA(:,7));
TIGHT = table2array(DATA(:,8));
NEWS = table2array(DATA(:,9));

% Defining the time interval of interest
Start = 62;
Stop = 182;

% Calculating the AR model variables of interest
Growth = (GDP(Start:Stop) - GDP(Start-4:Stop-4)) ./ GDP(Start-4:Stop-4) * 100;
Inflation = (CPI(Start:Stop) - CPI(Start-4:Stop-4)) ./ CPI(Start-4:Stop-4) * 100;
Rate = SIR(Start:Stop);
Spread = SPREAD(Start:Stop);
Vol = (VOL(Start:Stop) - VOL(Start-4:Stop-4)) ./ VOL(Start-4:Stop-4) * 100;
Variables = [Growth Inflation Rate Spread Vol];
K = size(Variables,2);

% Potential Instruments
Tight = TIGHT(Start:Stop);
News = NEWS(Start:Stop);

% Selecting the instrument
Instrument = Tight; % Textual measure of credit supply shock

% Select variable instrumented
inst = 4; % Credit spread

% Defining the lag
% [lag, AIC, BIC, HQC] = var_lag(Variables,8);
lag = 4; % selected on the basis of Mumtaz et al. (2018)

% VAR Estimation (OLS)
[BETA,SIGMA,ERR,X] = var_est(Variables,lag);

%% Identification of C matrix %%

% Identifying the relevant column of the "C" matrix (converts credit supply shock 
% into reduced form shocks)
Gamma = (ERR'*Instrument(lag+1:end))/length(Instrument(lag+1:end));

% Normalizing the relevant column of the "C" matrix (shock leads to 1 unit 
% increase of instrumented variable)
C_inst = Gamma / Gamma(inst,1);

%% Simple IRF %%
horizon = 40;

% Building the companion matrix
A = zeros(K*lag);
A(1:K,:) = BETA(2:end,:)';
A = A + diag(ones((lag-1)*K,1),-K);

% Calculation
IRF = zeros(horizon+1,K);
Shock = zeros(K*lag,1);
Shock(1:K,1) = C_inst;
for i = 1 : horizon+1
    response = Shock'*(A^(i-1))';
    IRF(i,:) = response(1,1:K);
end

figure
subplot(2,3,1)
hold on
plot(0:1:horizon,IRF(1:horizon+1,1))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("GDP Growth Response to Credit Shock")
hold off
subplot(2,3,2)
hold on
plot(0:1:horizon,IRF(1:horizon+1,2))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Inflation Response to Credit Shock")
hold off
subplot(2,3,3)
hold on
plot(0:1:horizon,IRF(1:horizon+1,3))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Policy Rate Response to Credit Shock")
hold off
subplot(2,3,4)
hold on
plot(0:1:horizon,IRF(1:horizon+1,4))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Rate Spread Response to Credit Shock")
hold off
subplot(2,3,5)
hold on
plot(0:1:horizon,IRF(1:horizon+1,5))
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Volume Growth Response to Credit Shock")
hold off

%% IRF with confidence intervals (bootstrap) %%

% Bootstrap
replications = 1000;
alpha = 0.1;
vector = 1:1:size(ERR,1);
reshuffled_index = RESHUFFLE(vector',10*replications);
ERR_RESHUFFLE = zeros(size(ERR,1),size(ERR,2),replications);
Instrument_RESHUFFLE = zeros(size(Instrument,1),replications);
j = 1;
k = 1;
while j <= 10*replications && k<= replications
    Instrument_RESHUFFLE_TENTATIVE(1:lag,1) = zeros(lag,1);
    for i = 1 : size(ERR,1)
        ERR_RESHUFFLE_TENTATIVE(i,:) = ERR(reshuffled_index(i,1,j),:);
        Instrument_RESHUFFLE_TENTATIVE(i+lag,1) = Instrument(reshuffled_index(i,1,j)+lag); 
    end
    if corr(ERR_RESHUFFLE_TENTATIVE(:,4),Instrument_RESHUFFLE_TENTATIVE(lag+1:end,1)) > 0.25
        ERR_RESHUFFLE(:,:,k) = ERR_RESHUFFLE_TENTATIVE;
        Instrument_RESHUFFLE(:,k) = Instrument_RESHUFFLE_TENTATIVE;
        k = k+1;
    end
    j = j+1;
end

for k = 1 : replications
    Y_star = zeros(size(Variables));
    Y_star(1:lag,:) = Variables(1:lag,:);
    for i = lag+1 : size(Variables,1)
        for j = 1 : lag
            Y_star(i,:) = Y_star(i,:) + Y_star(i-j,:) * BETA(2+(j-1)*K:1+j*K,:);
        end
        Y_star(i,:) = Y_star(i,:) + BETA(1,:) + ERR_RESHUFFLE(i-lag,:,k);
    end
    [BETA_star(:,:,k),SIGMA_star,ERR_star,X_star] = var_est(Y_star,lag);
    C_inst_star(:,k) = (ERR_RESHUFFLE(:,:,k)'*Instrument_RESHUFFLE(lag+1:end,k))/length(Instrument_RESHUFFLE(lag+1:end,k));
    C_inst_star(:,k) = C_inst_star(:,k) / C_inst_star(inst,k);
end

% IRF Calculation
horizon = 40;

IRF_star = zeros(horizon+1,K,replications);
for k = 1 : replications
    % Building the companion matrix
    A = zeros(K*lag);
    A(1:K,:) = BETA_star(2:end,:,k)';
    A = A + diag(ones((lag-1)*K,1),-K);
    % Initialization
    Shock = zeros(K*lag,1);
    Shock(1:K,1) = C_inst_star(:,k);
    % Loop on time steps
    for i = 1 : horizon+1
        response = Shock'*(A^(i-1))';
        IRF_star(i,:,k) = response(1,1:K);
    end
end

impact = 0.1 ; % Impact of initial credit shock on instrumented variable
 
IRF_new = impact * sort(IRF_star,3,'ascend');
IRF_low = IRF_new(:,:,replications * alpha/2);
IRF_med = IRF_new(:,:,replications * 0.5);
IRF_high = IRF_new(:,:,replications * (1-alpha/2));

figure 
subplot(2,3,1)
hold on
plotx2(0:horizon,[IRF_low(:,1) IRF_high(:,1)])
plot(0:1:horizon,IRF_med(:,1),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("GDP Growth Response to Credit Supply Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,2)
hold on
plotx2(0:horizon,[IRF_low(:,2) IRF_high(:,2)])
plot(0:1:horizon,IRF_med(:,2),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Inflation Response to Credit Supply Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,3)
hold on
plotx2(0:horizon,[IRF_low(:,3) IRF_high(:,3)])
plot(0:1:horizon,IRF_med(:,3),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Policy Rate Response to Credit Supply Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,4)
hold on
plotx2(0:horizon,[IRF_low(:,4) IRF_high(:,4)])
plot(0:1:horizon,IRF_med(:,4),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Rate Spread Response to Shock")
legend("Confidence interval (90%)","Median")
hold off
subplot(2,3,5)
hold on
plotx2(0:horizon,[IRF_low(:,5) IRF_high(:,5)])
plot(0:1:horizon,IRF_med(:,5),'color','red','LineWidth',2.0)
plot(0:1:horizon,zeros(horizon+1),'--','color','black')
title("Credit Volume Growth Response to Credit Shock")
legend("Confidence interval (90%)","Median")
hold off

%% Historical Decomposition %%

% Retrieving structural shock on instrumented variable
epsilon_inst = ERR * inv(SIGMA) * Gamma;

% Normalizing the structural shock on instrumented variable
epsilon_inst = epsilon_inst / sqrt(C_inst' * inv(SIGMA) * C_inst);

figure
plot(DATE(Stop-size(epsilon_inst)+1:Stop),epsilon_inst)
title("Structural Credit Shocks (normalized to epsilon/sigma)")

% Counterfactual if only credit supply shocks had occured
Var_CF = zeros(size(Variables,1),size(Variables,2));
for i = lag+1 : size(Var_CF,1)
    X = 0;
    for j = 1 : lag
        X = [X Var_CF(i-j,:)];
    end
    Var_CF(i,:) = X*BETA + ERR(i-lag,:)*inv(SIGMA)*Gamma*Gamma'/(Gamma'*inv(SIGMA_star)*Gamma); 
end

figure
subplot(2,1,1)
plot(DATE(Stop-size(Var_CF,1)+1:Stop),Var_CF(:,1),'--')
hold on
plot(DATE(Stop-size(Variables,1)+1:Stop),Variables(:,1))
plot(DATE(Stop-size(Variables,1)+1:Stop),zeros(length(Variables(:,1))),"--","color","black")
title("GDP Growth (YoY %)")
legend("Credit Supply Shock Induced","Historical")
hold off
subplot(2,1,2)
plot(DATE(Stop-size(Var_CF,1)+1:Stop),Var_CF(:,2),'--')
hold on
plot(DATE(Stop-size(Variables,1)+1:Stop),Variables(:,2))
plot(DATE(Stop-size(Variables,1)+1:Stop),zeros(length(Variables(:,2))),"--","color","black")
title("Inflation (YoY %)")
legend("Credit Supply Shock Induced","Historical")
hold off

%% Forecast Error Variance Decomposition %%

% Parameter setting
horizon = 40;

% Building the companion matrix
I = eye(lag*K);
A = [BETA(2:end,:)';I(1:end-K,:)];

% Phi matrices
Phi = zeros(K,K,horizon+1);
for i = 1 : horizon+1
    matrix = A^(i-1);
    Phi(:,:,i) = matrix(1:K,1:K);
end

% Calculating the MSE
MSE = zeros(K,K,horizon+1);
MSE(:,:,1) = SIGMA;
for i = 1 : horizon
    MSE(:,:,i+1) = MSE(:,:,i+1) + Phi(:,:,i+1) * MSE(:,:,1) * Phi(:,:,i+1)';
end

% Calculating the MSE due to structural credit shock only
MSE_cs = zeros(K,K,horizon+1);
MSE_cs(:,:,1) = Gamma*Gamma'/(Gamma'*inv(SIGMA)*Gamma);
for i = 1 : horizon
    MSE_cs(:,:,i+1) = MSE_cs(:,:,i+1) + Phi(:,:,i+1) * MSE_cs(:,:,1) * Phi(:,:,i+1)';
end

FEVD = zeros(horizon+1,K);
for i = 1 : K
    FEVD(:,i) = MSE_cs(i,i,:) ./ MSE(i,i,:)*100;
end

figure
subplot(2,3,1)
plot(0:horizon,FEVD(:,1))
title("FEVD : Credit Supply Shock on GDP Growth")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,2)
plot(0:horizon,FEVD(:,2))
title("FEVD : Credit Supply Shock on Inflation")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,3)
plot(0:horizon,FEVD(:,3))
title("FEVD : Credit Supply Shock on Policy Rate")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,4)
plot(0:horizon,FEVD(:,4))
title("FEVD : Credit Supply Shock on Credit Spread")
xlabel("Horizon")
ylabel("% MSE")
subplot(2,3,5)
plot(0:horizon,FEVD(:,5))
title("FEVD : Credit Supply Shock on Credit Volume Growth")
xlabel("Horizon")
ylabel("% MSE")

%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RESHUFFLED_SAMPLES] = RESHUFFLE(Vector,Number)
% This function reshuffles a vector (draw with replacement)

% Input :
% Vector : The vector to reshuffle
% Number : The number of new samples to produce

% Output :
% RESHUFFLED_SAMPLES : A matrix which columns are composed by new samples
% reshuffled from the original one.

[l,m] = size(Vector);
i = 1;
RESHUFFLED_SAMPLES = zeros(l,m,Number);
while i <= Number
    for k = 1:m
        for j=1:l
            RESHUFFLED_SAMPLES(j,k,i) = Vector(randi([1 l]),k);
        end
    end
    i = i+1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Beta,Sigma,ERR,x]=var_est(X,lag)

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
