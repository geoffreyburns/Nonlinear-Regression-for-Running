%% Nonlinear Regression Spring-Mass Analysis

% Author:
% Geoffrey Burns
% Michigan Performance Research Laboratory
% University of Michigan
% May 2020

% Files:
% vgrf_dat: vGRF data - 10 steps from subject C
% step_info: data to generate starting values
% 
% In-line Functions:
% PSvGRF: Parameterized Sinusoidal vGRF function
% PSvGRF_NLR: PSvGRF function for use with nonlinear regression tool
% PSvGRF_NLR_ME: PSvGRF function for use with mixed-effetcs nonlinear regression tool
% PSvGRF_NLR_3P: 3-parameter PSvGRF function for use with nonlinear regression tool
% PSvGRF_NLR_ME_3P: 3-parameter PSvGRF function for use with mixed-effetcs nonlinear regression tool
%
% Function Files:
% aL_estimate: Estimate atd and L0 from from non-vertical displacement parameter (A) and speed (see Appendix A)
%% Load Data
% vGRF data - vgrf_dat:
% [1,         2,       3    4,   5, ]
% [StepNumber StepTime vGRF Mass Leg]
%
% start value data - step_info:
% [1,         2, 3, 4, 5,   ]
% [StepNumber L0 tc tf Speed]
% 
% Note: Only the StepNumber, StepTime, vGRF, and Mass are necessary to run the NLR model. 
% The measured leg length (L0), contact time (tc), flight time (tf), and Speed are included here to calucate the model's initial seed values.
% This is not necessary, but it may be a good way of standardizing initial conditions so as not to bias the results.
% The Leg code (0 = L, 1 = R) is included here for reference

load 'NLR_demo_dat.mat'

%% Generate Starting Values
% As mentioned above, we would recommend standardizing the starting values for the NLR
% optimizer with tradional spring-mass calculations (with L0 either
% measured or estimated from height; contact time observed; touchdown angle estimated from leg
% length, tc, and average speed; and stiffness from either McMahon and
% Cheng (1990) or Morin et al. (2005)). 
% The example below is a calculation via Morin et al. (2005):

startvals_L0 = step_info(:,2);
startvals_t0 = step_info(:,3);
startvals_a0 = acos((step_info(:,5).*step_info(:,3)/2)./step_info(:,2));
startvals_k0 = (vgrf_data(1,4).*9.80665.*(pi/2).*(step_info(:,4)./step_info(:,3)+1))...
    ./((step_info(:,2)-sqrt(step_info(:,2).^2-(step_info(:,5).*step_info(:,3)./2).^2))...
    + ((vgrf_data(1,4).*9.80665.*(pi/2).*(step_info(:,4)./step_info(:,3)+1)).*step_info(:,3).^2)./(vgrf_data(1,4).*pi^2)-9.80665/8.*(step_info(:,3)).^2);
startvals = [startvals_k0 startvals_a0 startvals_L0 startvals_t0];

%% Nonlinear Regression Function

PSvGRF_NLR = @(b,X)((b(1).*(b(3)-(9.80665.*(b(4)).^2)./8-b(3).*sin(b(2))))...
    ./(1-(b(1).*(b(4)).^2)./(X(:,2).*pi^2)))...
    .*sin(X(:,1).*(pi./(b(4))))...
    .*(1+(-1./(1+exp(-1e10.*(X(:,1)-(b(4)))))));

%% Nonlinear Regression Model

time    = [vgrf_data(:,2), vgrf_data(:,4)];
force   = vgrf_data(:,3);
beta0 = mean(startvals);

sm_nlr = fitnlm(time, force, PSvGRF_NLR, beta0)

% Note: Model details and objects (Residuals, Diagnostics, Error Terms,
% Fitting Criterion, etc.) can be explored in the "sm_nlr" object output

%% Plot the NLR Model

beta_nlr = table2array(sm_nlr.Coefficients(1:4,1))';
step    = vgrf_data(:,1);

clf

hp = gscatter(time(:,1),force,step,[],[],10);
h.MarkerFaceAlpha = .4;
xlabel('Time (secs)')
ylabel('Force (N)')
title('{\bf Spring-Mass vGRF via NLR}')
legend('off')
xx = linspace(0,.2)';
xlim([0 .21]);
mass = time(1,2);
line(xx, PSvGRF(beta_nlr, xx, mass),'linewidth',2,'color','k');


%% Set Constraints
% Note: the nonlinear optimizers used by MatLab do not allow for
% constraints. To remedy this, we added logistic multipliers for each
% parameter on desired upper and lower bounds so that the value is 
% multiplied by 1 within the bounds and 1e30 outside the bounds, 
% effectively constraining the objective function. The constraints used
% below are quite broad from a physiological perspective (an individual is
% very unlikely to be able to comfrotably run at the given speed with any
% of these upper or lower bounds), and the researcher can use tighter
% ranges if he or she deems it warranted.

k_lb = 5000;
k_ub = 20000;
a_lb = 1.1;
a_ub = 1.3;
l_lb = 0.80;
l_ub = 1.20;
t_lb = 0.12;
t_ub = 0.40;

%% Mixed-Effects Nonlinear Regression Function

PSvGRF_NLR_ME = @(B,X)((B(:,1).*(B(:,3)-(9.80665.*(B(:,4)).^2)./8-B(:,3).*sin(B(:,2))))...
    ./(1-(B(:,1).*(B(:,4)).^2)./(X(:,2).*pi^2)))...
    .*sin(X(:,1).*(pi./(B(:,4))))...
    .*(1+(-1./(1+exp(-1e10.*(X(:,1)-(B(:,4)))))))...
    +(1/(1+exp(-1e10.*(B(:,1)-k_ub)))+1/(1+exp(-1e10.*(k_lb-B(:,1))))).*1e10...
    +(1/(1+exp(-1e10.*(B(:,2)-a_ub)))+1/(1+exp(-1e10.*(a_lb-B(:,2))))).*1e10...
    +(1/(1+exp(-1e10.*(B(:,3)-l_ub)))+1/(1+exp(-1e10.*(l_lb-B(:,3))))).*1e10...
    +(1/(1+exp(-1e10.*(B(:,4)-t_ub)))+1/(1+exp(-1e10.*(t_lb-B(:,4))))).*1e10;

%% Set NLR Options and ME Starting Values

% Inital Parameter Set

beta0 = mean(startvals);

% Initial Covariance Matrix for Mixed-effects
% Note: Here, we used a conservative estimate of twice the variance of each
% of the starting parameters on a diagonal identity matrix (and the variance 
% of the leg length was estimated via the touchdown angle variance with the 
% mean contact time). The structure of this can be maniuplated depending on 
% the research question and assumptions about the relationships of the
% paramters

cov0 = (2*var(startvals)).*eye(size(beta0,2));
cov0(3,3) = 2*var(.5*step_info(1,2)*beta0(4)./cos(startvals(:,2)));

options = statset('Display', 'iter', 'FunValCheck', 'off', ...
                    'DerivStep', [1 1e-5 1e-5 1e-5], 'TolFun', 1e-2, ...
                    'OutputFcn', @nlmefitoutputfcn);

time    = [vgrf_data(:,2), vgrf_data(:,4)];
force   = vgrf_data(:,3);
step    = vgrf_data(:,1);

%% Performing the Nonlinear Regression

[beta_nlme, PSI, stats_nlme, B] = nlmefitsa(time, force, step, ...
                            [], PSvGRF_NLR_ME, beta0, ...
                            'NIterations', [300 200 100], ...
                            'NBurnIn', 5, 'NChains', 5, 'Options', options, ...
                            'Cov0', cov0, ...
                            'REParamsSelect', [1 2 3 4], ...
                            'ComputeStdErrors', true, 'LogLikMethod', 'lin',...
                            'OptimFun', 'fminsearch');
            
%% Plotting - NLR ME
%%
% Fixed Effect

clf
hp = gscatter(time(:,1),force,step,[],[],10);
h.MarkerFaceAlpha = .4;
xlabel('Time (secs)')
ylabel('Force (N)')
title('{\bf Spring-Mass vGRF via ME NLR}')
legend('off')
xx = linspace(0,.2)';
xlim([0 .21]);
mass = time(1,2);
line(xx, PSvGRF(beta_nlme, xx, mass),'linewidth',2,'color','k');

%%
% Each Step's Random Effect
for j=1:10
    phir = (beta_nlme+B(:,j))';
    hl = line(xx,PSvGRF(phir, xx, mass),'linewidth',2,'color',hp(j).Color);
    hl.Color(4) = 0.6;
end

line(xx,PSvGRF(beta_nlme, xx, mass),'linewidth',3,'color','k');


%% Three Parameter Model 
% (Appendix A Alternative)

%% Three Parameter: Starting Values
% Generate new starting value for 3-parameter model
% A = L0 - L0*sin(a)

startvals_3P = [startvals(:,1) startvals(:,3)-startvals(:,3).*sin(startvals(:,2)) startvals(:,4)];

%% Three Parameter: Nonlinear Regression Function

PSvGRF_NLR_3P = @(b,X)((b(1).*(b(2)-(9.80665.*(b(3)).^2)./8))...
    ./(1-(b(1).*(b(3)).^2)./(X(:,2).*pi^2)))...
    .*sin(X(:,1).*(pi./(b(3))))...
    .*(1+(-1./(1+exp(-1e10.*(X(:,1)-(b(3)))))));

%% Three Parameter: Nonlinear Regression Model

time    = [vgrf_data(:,2), vgrf_data(:,4)];
force   = vgrf_data(:,3);
beta0 = mean(startvals_3P);

sm_nlr_3P = fitnlm(time, force, PSvGRF_NLR_3P, beta0)

beta_nlr_3P = table2array(sm_nlr_3P.Coefficients(1:3,1))';


% Estimate L0 and atd from A and speed (Appendix A)
[~,~,~,~,~,L0_est, atd_est] = aL_estimate(beta_nlr_3P(1), beta_nlr_3P(2), beta_nlr_3P(3), mass, step_info(1,5))



%% Three Parameter: ME NLR Starting Values

A_lb = l_lb-l_lb*sin(a_ub);
A_ub = l_ub-l_ub*sin(a_lb);

beta0_3P = mean(startvals_3P);
cov0_3P = (2*var(startvals_3P)).*eye(size(beta0_3P,2));


%% Three-Parameter: Mixed-Effects Nonlinear Regression Function

% ME NLR Function
PSvGRF_NLR_ME_3P = @(B,X)((B(:,1).*(B(:,2)-(9.80665.*(B(:,3)).^2)./8))...
    ./(1-(B(:,1).*(B(:,3)).^2)./(X(:,2).*pi^2)))...
    .*sin(X(:,1).*(pi./(B(:,3))))...
    .*(1+(-1./(1+exp(-1e10.*(X(:,1)-(B(:,3)))))))...
    +(1/(1+exp(-1e20.*(B(:,1)-k_ub)))+1/(1+exp(-1e10.*(k_lb-B(:,1))))).*1e30...
    +(1/(1+exp(-1e20.*(B(:,2)-A_ub)))+1/(1+exp(-1e20.*(A_lb-B(:,2))))).*1e30...
    +(1/(1+exp(-1e20.*(B(:,3)-t_ub)))+1/(1+exp(-1e20.*(t_lb-B(:,3))))).*1e30;

% ME NLR Model
[beta_nlme_3P, PSI_3P, stats_nlme_3P, B_3P] = nlmefitsa(time, force, step, ...
                            [], PSvGRF_NLR_ME_3P, beta0_3P, ...
                            'NIterations', [300 200 100], ...
                            'NBurnIn', 5, 'NChains', 5, 'Options', options, ...
                            'Cov0', cov0_3P, ...
                            'REParamsSelect', [1 2 3], ...
                            'ComputeStdErrors', true, 'LogLikMethod', 'lin',...
                            'OptimFun', 'fminsearch');

% Estimate L0 and a_td for each step from Random Effects (Appendix A)

phir_3P_est = zeros(10,4);

for j=1:10
    phir_3P(j,:) = (beta_nlme_3P+B_3P(:,j))';
end

 for j = 1:10
    [~,~,~,~,~,L0_est, atd_est] = aL_estimate(phir_3P(j,1), phir_3P(j,2), phir_3P(j,3), mass, step_info(1,5));
    
    phir_3P_est(j,:) = [phir_3P(j,1), L0_est, atd_est, phir_3P(j,3)];
 end


%% PS vGRF Function

function PSvGRF = PSvGRF(b, t, m)
    
k = b(1);
a = b(2);
L0 = b(3);
CT = b(4);

PSvGRF =    k.*(L0-(9.80665*CT^2)/8-L0.*sin(a))./(1-(k.*CT^2)./(m.*pi^2)) ...
                .*sin(t(:,1).*(pi./CT)) ...
                .*(1+(-1./(1+exp(-1e10.*(t(:,1)-CT))))); 
end

