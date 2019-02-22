% partA task a

minSize = min([numel(SPY) numel(GOVT) numel(EEMV)]);
SPY = SPY(numel(SPY)-minSize+1:end);
GOVT = GOVT(numel(GOVT)-minSize+1:end);
EEMV = EEMV(numel(EEMV)-minSize+1:end);

% return of asset
r_SPY=diff(SPY)./SPY(1:end-1,:);
r_GOVT=diff(GOVT)./GOVT(1:end-1,:);
r_EEMV=diff(EEMV)./EEMV(1:end-1,:);

% expected return
OnePlusSPY = 1+ r_SPY;
ProdOnePlusSPY = prod(OnePlusSPY);
u_SPY = power(ProdOnePlusSPY, (1/numel(OnePlusSPY))) - 1;

OnePlusGOVT = 1+ r_GOVT;
ProdOnePlusGOVT = prod(OnePlusGOVT);
u_GOVT = power(ProdOnePlusGOVT, (1/numel(OnePlusGOVT))) - 1;


OnePlusEEMV = 1+ r_EEMV;
ProdOnePlusEEMV = prod(OnePlusEEMV);
u_EEMV = power(ProdOnePlusEEMV, (1/numel(OnePlusEEMV))) - 1;

% standard deviation
stdSPY = std(r_SPY);
stdGOVT = std(r_GOVT);
stdEEMV = std(r_EEMV);

% covariance
cov_SPY_GOVT_m = cov(r_SPY, r_GOVT);
cov_SPY_GOVT=cov_SPY_GOVT_m(1,2);
cov_SPY_EEMV_m = cov(r_SPY, r_EEMV);
cov_SPY_EEMV=cov_SPY_EEMV_m(1,2);
cov_EEMV_GOVT_m = cov(r_EEMV, r_GOVT);
cov_EEMV_GOVT=cov_EEMV_GOVT_m(1,2);

%% 

% task b Generate an efficient frontier between three assets

A = -[u_SPY, u_EEMV, u_GOVT];
Q = [stdSPY^2, cov_SPY_EEMV, cov_SPY_GOVT;
    cov_SPY_EEMV, stdEEMV^2, cov_EEMV_GOVT;
    cov_SPY_GOVT, cov_EEMV_GOVT, stdGOVT^2;];
c = [0 0 0]';
Aeq=[1 1 1]; 
beq=[1];
lb=[0;0;0;]; % without short selling
ub=[inf;inf;inf;];
lb_s=[-inf;-inf;-inf];

R=[0];
i = 0;
total_dot = 20;
R_series = zeros(total_dot,1);
rho_series = zeros(total_dot,1);
rho_series_s = zeros(total_dot,1);
weight_series = zeros(total_dot,3);
weight_series_s = zeros(total_dot,3);
while i < total_dot
    i = i + 1;
    R_series(i,1)= R;
    b = -R;
    
    [x, fval] = quadprog(Q,c,A,b,Aeq,beq,lb,ub);
    rho_series(i,1) = sqrt(fval);
    weight_series(i,:) = x;
    
    [x, fval] = quadprog(Q,c,A,b,Aeq,beq,lb_s,ub);
    rho_series_s(i,1) = sqrt(fval);
    weight_series_s(i,:) = x;
    
    R = R + 0.0005;
    
end 

%Plotting Graph
figure; hold on
grid on
grid minor
a1 = plot(rho_series,R_series,'-*');
a2 = plot(rho_series_s,R_series,'-o');
title('The efficient frontier of MVO');
xlabel('Volatility');
ylabel('Expected return goal R');
legend([a1;a2], "EF without short selling", "EF with short selling")
hold off
%% 
% part B 

% the order of stocks in this session is 
% 1. SPY, 2. GOVT, 3. EEMV, 4. ACN, 5. BR, 6. CBOE, 7. ICE, 8. CME

% ACN = table2array(ACN);
% BR = table2array(BR);
% CBOE = table2array(CBOE);
% ICE = table2array(ICE);
% CME = table2array(CME);

minSize = min([numel(SPY) numel(GOVT) numel(EEMV) numel(ACN) numel(BR) numel(CBOE) numel(ICE) numel(CME)]);
raw_series = zeros(minSize,8);
raw_series(:,1) = SPY(numel(SPY)-minSize+1:end);
raw_series(:,2) = GOVT(numel(GOVT)-minSize+1:end);
raw_series(:,3) = EEMV(numel(EEMV)-minSize+1:end);
raw_series(:,4) = ACN(numel(ACN)-minSize+1:end);
raw_series(:,5) = BR(numel(BR)-minSize+1:end);
raw_series(:,6) = CBOE(numel(CBOE)-minSize+1:end);
raw_series(:,7) = ICE(numel(ICE)-minSize+1:end);
raw_series(:,8) = CME(numel(CME)-minSize+1:end);

% return of asset
r_series = zeros(minSize-1,8);
for i=1:8
    currnet_series = raw_series(:,i);
    r_series(:,i) = diff(currnet_series)./currnet_series(1:end-1,:);
end

% expected return
u_series = zeros(1,8);
for i=1:8
    currnet_series = r_series(:,i);
    one_plus = 1+ currnet_series;
    prod_one_plus = prod(one_plus);
    u_series(1,i) = power(prod_one_plus, (1/numel(one_plus))) - 1;
end

% standard deviation
std_series = std(r_series);

% covariance
cov_series = zeros(8,8);
for i=1:8
    for j=1:8
    cov_temp = cov(r_series(:,i), r_series(:,j));
    cov_series(i,j)=cov_temp(1,2);
    end
end

A = -u_series;
Q = cov_series;
c = [0 0 0 0 0 0 0 0]';
Aeq=[1 1 1 1 1 1 1 1]; 
beq=[1];
lb=[0;0;0;0;0;0;0;0;]; % without short selling
ub=[inf;inf;inf;inf;inf;inf;inf;inf;];
lb_s=[-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;];

R=[0];
i = 0;
total_dot = 20;
R_series = zeros(total_dot,1);
rho_series_2 = zeros(total_dot,1);
rho_series_s_2 = zeros(total_dot,1);
weight_series = zeros(total_dot,8);
weight_series_s = zeros(total_dot,8);
while i < total_dot
    i = i + 1;
    R_series(i,1)= R;
    b = -R;
    
    [x, fval] = quadprog(Q,c,A,b,Aeq,beq,lb,ub);
    rho_series_2(i,1) = sqrt(fval);
    weight_series(i,:) = x;
    
    [x, fval] = quadprog(Q,c,A,b,Aeq,beq,lb_s,ub);
    rho_series_s_2(i,1) = sqrt(fval);
    weight_series_s(i,:) = x;
    
    R = R + 0.0005;
    
end 

%Plotting Graph
figure; hold on
grid on
grid minor
a3 = plot(rho_series_2,R_series,'-*');
a4 = plot(rho_series_s_2,R_series,'-o');
title('The efficient frontier of MVO');
xlabel('Volatility');
ylabel('Expected return goal R');
legend([a1;a2], "EF without short selling", "EF with short selling")
hold off

figure; hold on
grid on
grid minor
a1 = plot(rho_series,R_series,'-*');
a2 = plot(rho_series_s,R_series,'-o');
a3 = plot(rho_series_2,R_series,'-*');
a4 = plot(rho_series_s_2,R_series,'-o');
title('The efficient frontier of MVO comparison');
xlabel('Volatility');
ylabel('Expected return goal R');
legend([a1;a2;a3;a4], "EF for the first three stock without short selling", "EF for the first three stock with short selling","EF for all eight stock without short selling", "EF for all eight stock with short selling");
hold off
