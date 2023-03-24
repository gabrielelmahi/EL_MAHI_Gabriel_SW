% Advanced Macroeconomics I - Assignment 3
% Group 10
% Credit Market Shock and Real Activity
% Method 3: Sign restrictions Approach

close all;
clear;
clc;

%inputs
data_or = readtable('Data.xlsx','Sheet','DATA');

Spread = table2array(data_or(:,4));
Policy_Rate = table2array(data_or(:,3));
CPI = table2array(data_or(:,5));
Real_GDP = table2array(data_or(:,6));
Credit_Supply = table2array(data_or(:,7));


GDP_Growth = (Real_GDP(5:end) - Real_GDP(1:end-4)) ./ Real_GDP(1:end-4) * 100;
Inflation = (CPI(5:end) - CPI(1:end-4)) ./ CPI(1:end-4) * 100;
Policy_Rate = Policy_Rate(5:end);
Spread = Spread(5:end);
Credit_Supply_Growth = (Credit_Supply(5:end) - Credit_Supply(1:end-4)) ./ Credit_Supply(1:end-4) * 100;

% the data matrix
data=[GDP_Growth Inflation Policy_Rate Spread Credit_Supply_Growth];

[lag, AIC, BIC, HQC] = var_lag(data,10);

reps = 1000; % Total number of replications
update = 1; % If update = 1, then print every update iteration
horizon = 40; % The Forecast horizon
Lags = 4;   % number of lags in the VAR


%Pattern : 0=no restriction, 1 means >0 and -1 <0
pattern=[0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0
    -1 -1 -1 1 -1];

% timer start
tic

Y = data;
N=size(Y,2);
%take lags
X=[];
for j=1:Lags
    X=[X lag0(data,j) ];
end

X=[X ones(size(X,1),1)];

Y=Y(Lags+1:end,:);
X=X(Lags+1:end,:);

T=size(Y,1);

fsave=zeros(reps,N,horizon,N);

% Display empty line (to separate samples in the screen shot)
disp(sprintf(' '))


%% Estimation
% OLS estimation of reduced-form coefficients
beta0 = inv(X'*X)*X'*Y;
RESIDUALS = Y-X*beta0;
sigma = RESIDUALS'*RESIDUALS/(T-Lags);
Ctilde=chol(sigma);

i=1;
j=1;

% Tant que le nombre de réplication n'est pas atteint
while j<reps+1

    % Display progress:
    if update==1
        disp(sprintf(' Replication %s of %s. / Total draws %s', ...
            num2str(j), num2str(reps),num2str(i)) );
    end

    % chck is negative as long as we don't find a valid candidate matrix
    chck=-1;

    while chck<0
        % i : number of total draws
        i = i + 1;
        
        % Proceed to QR decomposition of a random matrix K which the elements
        % are drawn from a standard normal distribution
        K = normrnd(0,1,N,N);
        [Q,R] = qr(K);

        C=(Q*Ctilde);  %candidate draw

        % Finding a matrix respecting sign restrictions constraints
        C_New=geta0(C,pattern,beta0,Lags);

        % If the candidate matrix is valid. i.e non null
        if sum(sum(C_New))~=0
            chck = 1;

            irfmat=zeros(N,horizon,N);

            for jj=1:N
                shock=zeros(1,N);
                shock(jj)=1;
                irfmat(jj,:,:)=impulse(beta0,N,Lags,C_New,shock,horizon+Lags);
            end

            fsave(j,:,:,:)=irfmat;

            % Proceed to forecast error descomposition for each replication
            fevd= FEVD(horizon,N,Lags,beta0,C_New);
            fevdmat(j,:,:)=fevd;


            % Proceed the historical decomposition for each replication
            Aols=[beta0(end,:);beta0(1:end-1,:)]';
            errormat = RESIDUALS(1:size(Y,1),:)*inv(C_New);

            histmat(j,:,:,:)=Historical_Decomposition(Aols,C_New',data,...
                errormat');

            % increment the number of successful draws
            j=j+1;
        end
    end
end

% timer stop
toc

%plot the response to shock 3 the policy shock

confidence_interval = [6 94];
%plot the response to shock 5 the credit supply shock
temp=squeeze(fsave(:,5,:,:));
tt=0:horizon-1;
xx = zeros(1,horizon);


%% Plot IRFs

% NORMALISATION 
% Median response of the spread reesponse to credit supply shock
mediane_4 = median(temp(:,:,4));
normalisation = mediane_4(1)/0.1;

varNames = {'GDP growth','Inflation','Policy Rate','Spread','Credit volume growth'};

    figure()
    for k=1:N
        subplot(2,3,k)
        temp_conf_int=squeeze(prctile(temp(:,:,k)/normalisation,confidence_interval,1));
        plotx2(tt,temp_conf_int')
        hold on
        mediane_1 = median(temp(:,:,k)/normalisation);
        plot(tt,mediane_1,'--','color','blue','LineWidth',2)
        hold on
        plot(tt,xx,'color','black','LineWidth',1)
        title(['Credit shock on ' varNames{k}])
        legend('Confidence interval - 90%','Median','Location','northeast');
    end

    
    figure()
    temp_conf_int=squeeze(prctile(temp(:,:,1)/normalisation,confidence_interval,1));
    plotx2(tt,temp_conf_int')
    hold on
    mediane_1 = median(temp(:,:,1)/normalisation);
    plot(tt,mediane_1,'--','color','blue','LineWidth',2)
    hold on
    plot(tt,xx,'color','black','LineWidth',1)
    legend('Confidence interval - 90%','Median','Location','southeast');
    title(['Credit shock on ' varNames{1}])

    figure()
    temp_conf_int=squeeze(prctile(temp(:,:,2)/normalisation,confidence_interval,1));
    plotx2(tt,temp_conf_int')
    hold on
    mediane_1 = median(temp(:,:,2)/normalisation);
    plot(tt,mediane_1,'--','color','blue','LineWidth',2)
    hold on
    plot(tt,xx,'color','black','LineWidth',1)
    legend('Confidence interval - 90%','Median','Location','southeast');
    title(['Credit shock on ' varNames{2}])


%% Plot Historical decompositions

% Plot historical decomposition of GDP growth and Inflation with regards to
% credit supply shocks

% Evolution of inflation only with credit supply shocks - Counterfactual
GDP_Growth_5=squeeze(prctile(histmat(:,:,5,1),[50 6 94],1)); % 5 : Credit supply shock, 1 : impact on GDP growth

% Evolution of inflation only with credit supply shocks - Counterfactual
Inflation_5=squeeze(prctile(histmat(:,:,5,2),[50 6 94],1)); % 5 : Credit supply shock, 2 : impact on inflation

% define the considered period
period=1977;
freq = 4;
gend = size(GDP_Growth_5(1,:),2);
x = period-2/freq:(1/freq):period+(gend-2)/freq;

% Plot the historical decomposition of GDP Growth
figure()
plot(x(2:end),GDP_Growth(4:end),'LineWidth',2);
hold on
plot(x(2:end),GDP_Growth_5(1,:)','--','color','red','LineWidth',2);
hold on
plot(x(2:end),zeros(1,gend),'k-')
legend('Historical GDP growth','Contribution of credit supply shock','Location','northeast');
title('Historical decomposition of GDP growth')

% Plot the historical decomposition of Inflation
figure()
plot(x(2:end),Inflation(4:end),'LineWidth',2);
hold on
plot(x(2:end),Inflation_5(1,:)','--','color','red','LineWidth',2);
hold on
plot(x(2:end),zeros(1,gend),'k-')
legend('Historical Inflation','Contribution of credit supply shock','Location','northeast');
title('Historical decomposition of Inflation')


%% Plot FEVD

% Plot FEVD of GDP growth
tmp=squeeze(mean(fevdmat)); %the median may not add up to 1

figure()
% plot([tmp(:,1) tmp(:,6) tmp(:,11) tmp(:,16) tmp(:,21)])
plot(tmp(:,1),'color','red')
hold on
plot(tmp(:,6),'color','blue')
hold on
plot(tmp(:,11),'color','black')
hold on
plot(tmp(:,16),'color','magenta')
hold on
plot(tmp(:,21),'color','green','LineWidth',2)

legend('Shock 1 : GDP shock','shock2 : Inflation shock','shock3 : Policy rate shock','shock4 : Spread shock','shock5 : Credit supply shock','Location','best')
title('Contribution to Foreceast Error Variance - GDP Growth')


figure()
plot(tmp(:,2),'color','red')
hold on
plot(tmp(:,7),'color','blue')
hold on
plot(tmp(:,12),'color','black')
hold on
plot(tmp(:,17),'color','magenta')
hold on
plot(tmp(:,22),'color','green','LineWidth',2)

legend('Shock 1 : GDP shock','shock2 : Inflation shock','shock3 : Policy rate shock','shock4 : Spread shock','shock5 : Credit supply shock','Location','best')
title('Contribution to Foreceast Error Variance - Inflation')


%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

end

function y=impulse(b,n,l,v,s,t)
%b= VAR coefs
%n=number of variables
%l=lag length
%v=A0 matrix
%s=shock vector
%t=horizon

e=zeros(t+l,n);

e(l+1,:)=s;

y=zeros(t+l,n);

for k=l+1:t
    x=[];
    for i=1:l
        for j=1:n
            x=[x y(k-i,j)];
        end
    end
    y(k,:)=([x 0]*reshape(b,n*l+1,n))+(e(k,:)*v);
end
y=y(l+1:size(y,1)-l,:);
end


function a0new=geta0(a0old,pattern,beta,L)


n=size(a0old,1);
out=[];

% vector of dimension n
% equal o 1 if there is a sign restriction of shock n on variables
tempid= sum(abs(pattern),2)>0;

% Extract the non null pattern
patterni=pattern(tempid,:);

% Define the number of shocks to identify
nshocks=size(patterni,1);

% For each shock
for MM=1:nshocks

    patternj=patterni(MM,:); % Signs to check

    zmp=-1;

    % for each variable
    for j=1:n
        
        z1=a0old(j,:);

        % Compute the impulse response following a shock from variable j
        shock=zeros(1,n);
        shock(j)=1;

        % the response of variables to shock(j) in the first period
        irf=impulse(beta,n,L,a0old,shock,1+L);
        
        % Check if the irf respects the sign restrictions in the first
        % period
        check=check_restriction(irf,patternj);

        if check==1
            zmp=1*z1;
        end

    end

    % if the candidate matrix doesn't respect the sign restriction
    if zmp==-1
        out=[out;nan(1,n)];
    else % if the candidate matrix respects the sign restriction
        out=[out;zmp];
    end

end

% if for a shock the candidate matrix doesn't respect the sign restriction
% then we give a null matrix as an output
if sum(sum(isnan(out)))>0
    a0new=zeros(n,n);
else
    % out
    outx=[];

    % we keep a draw if we find exactly one corresponding shock
    for i=1:n
 	   zA=a0old(i,:);
       eam=[];

       for j=1:size(out,1)
           tempcheck= sum(zA == out(j,:));
    	   eam=[eam tempcheck];
       end
       eALL=sum(eam,2);
       if eALL==0
  	     outx=[outx;zA];
       end
    end

    %find appropriate row to insert shocks
    a0new=zeros(n,n);
    i=1;
    j=1;
    jj=1;
    for i=1:n
        tempcheck=sum(abs(pattern(i,:)),2)==0;
        if tempcheck==1
            a0new(i,1:n)=outx(j,:);
            j=j+1;
        else
            a0new(i,1:n)=out(jj,:);
            jj=jj+1;
        end
    end

end
end

function FEVDOUT = FEVD(HORIZON,N,L,betam,A0NEW)
% This function computes the contribution of each structural shock 
% to the forecast error decomposition of model variables

% Inputs :
%   - Horizon : number of periods
%   - N : Number of variables
%   - M : Number of lags
%   - betam : OLS estimates of beta
%   - A0NEW : Candidate matrix
% Outputs :
%   - FEVDOUT : matrix specifing the contribution of each structural shock
%               to fevd.

% Initialisation
FEVDOUT=0;


%forecast error descomposition
yhatall=zeros(HORIZON,N*N);
totvar=0;
jjj=1;
for iii=1:N
    yhat=zeros(HORIZON+L,N);
    vhat=zeros(HORIZON+L,N);
    vhat(L+1,iii)=1;

    for j=L+1:HORIZON+L
        xhat=[];
        for jj=1:L
            xhat=[xhat yhat(j-jj,:)];
        end
        xhat=[xhat 0];
        yhat(j,:)=xhat*reshape(betam,N*L+1,N)+vhat(j,:)*A0NEW;
    end
    yhatall(:,jjj:jjj+N-1)=yhat(L+1:end,:);
    totvar=totvar+cumsum(yhat(L+1:end,:).^2);
    jjj=jjj+N;
end

FEVDOUT=cumsum(yhatall.^2)./repmat(totvar,1,N);

end

function check=check_restriction(x,e)

% Inputs
%   - x : IRF : the dimension of x is equal to horizon*number of variables
%   - e : the pattern to respect. vector of dimension 1*number of variables
%   - timei : the time restriction, same as e.

% Output
%   - check = 1 means that the IRF respects the sign restriction. 0
%   otherwise


out=[];
% Number of variables 
n=size(x,2);

% Apply this loop to each the response of each variable 
for i=1:n
    % the impulse response of variable i to a shock
    x1=x(:,i);
    % the pattern corresponding to the reponse of variables i to a shock
    x2=e(i);
    
    out1=0;

    % If the impulse response of variable i has to respect a sign
    % restriction for a period x3
    if x2 ~= 0
        

        x1neg=sum(x1(1)<0) == 1;
        x1pos=sum(x1(1)>0) == 1;

        % if the response of variable i is negative and the pattern is also
        % negative
        if x1neg && x2<0
            out1=1;
        end
        % if the response of variable i is positive and the pattern is also
        % positive
        if x1pos && x2>0
            out1=1;
        end
    end
    out=[out out1];
end


% if the sum of out is equal to the number of restrictions in the pattern
% then 
check=sum(out,2)==sum(abs(e),2);

end


function yfs = Historical_Decomposition(Aols,IdentMat,datamat,...
    errormat)
% this function realizes the historical decomposition of each endogenous
% variable in the light of our model
% Inputs :
%   - Aols : our beta coefficients from the OLS estimation
%   - IdentMat : transpose of our Identified matrix
%   - datamat : matrix of data
%   - errormat : matrix of errors
% Output :
%   - yfs : the matrix of historical decomposition.

% dy : number of variables
% pdy : dimension of beta (dy*Lags) without the constant
[dy,pdy] = size(Aols(:,2:end));

% number of lags
Lags = pdy/dy;

% x = Y
% lx = X
[x,lx] = makelags(datamat',Lags);

% T : number of observations : time dimension
T = size(x,1);

% Compute the companion matrix of Aols_wc (betas without the constant)
Aols_wc = Aols(:,2:end);
[n_Aols_wc k_Aols_wc] = size(Aols_wc);
porder = k_Aols_wc/n_Aols_wc;
Acomp = [Aols_wc; [kron(eye(porder-1),eye(n_Aols_wc)) zeros((porder-1)*n_Aols_wc,n_Aols_wc)]];

% Realize the historical decomposition

TermD = zeros(pdy,dy);
TermD(1:dy,:) = IdentMat;
Jmat = [eye(dy) zeros(dy,(Lags-1)*dy)];

xfd = zeros(Lags*dy,T+1,dy);
yfd = zeros(dy,T+1,dy);

for j = 1 : dy
    termA  = zeros(dy,T+1);
    termA(j,2:end) = errormat(j,:);
    for i = 2 : T+1
        xfd(:,i,j) = Acomp * xfd(:,i-1,j)+ TermD * termA(:,i);
        yfd(:,i,j) = Jmat * xfd(:,i,j);
    end
end

yfs = zeros(T+1,dy,dy);
for i = 1 : dy
    for j = 1 : dy
        yfs(:,j,i)=yfd(i,:,j)';
    end
end
end

function out=lag0(x,p)
% this function lags the variable x and places zeros in the rows
% corresponding to the first p periods.
% inputs : 
%   - x : vector or matrix x = x=(a_{1}; a_{2}; ...; a_{n}) where a_{i}
%         are row vectors
% output : 
%   - out : (0; ... ; 0 ; a_{1}; a_{2}; ....; a_{n-p})

% Compute the number of rows and columns of input x
[R,C]=size(x);

% Take the first R-p rows of matrix x
x1=x(1:(R-p),:);

% Preceed them with p rows of zeros and return
out=[zeros(p,C); x1];

end


function [ydata, laggedydata] = makelags(ydata0,nL)

[Ty Ny]=size(ydata0);

if Ny > Ty
    ydata0 = ydata0';
    [Ty Ny]=size(ydata0);
end

ydata = ydata0(nL+1:Ty,:);

for i = 1 : nL
    laggedydata(:,(i-1)*Ny+1:i*Ny) = ydata0(nL+1-i:Ty-i,:);
end
end





