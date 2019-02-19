% Simulate counter-factual: If No Cannibalization
delete 'Diary_cf_NoCannibn'; diary('Diary_cf_NoCannibn');

global beta delta State Exit Adopt T Pi V EV Policy;
load CheckPi_NoCannibn.mat;

% Basic parameters
scaling = 1000;     % Re-scale profits, to avoid explosion due to exp(pi) or exp(EV)
beta = .80;         % Discount factor
delta = 1.10;      % Decline in SC (= MC_new's drop rate)
potent = 4;         % Number of potential entrants (assumption) for all years
State(:,4) = potent*ones(18,1);
T = size(State,1);  % Number of periods (years)

% Simulation parameters (MLE of beta = .80, npe =4, delta = 1.10)
phi = -0.147362421730779;
kappa_inc = 1.24394591815957;
kappa_ent = 2.25381667717391;
Theta = [phi, kappa_inc, kappa_ent];

fprintf('\n  -----------------------------------------');
fprintf('\n | STRATEGIC (RATIONAL-EXPECTATIONS) MODEL |');
fprintf('\n |         WITHOUT  CANNIBALIZATION        |');
fprintf('\n  -----------------------------------------');

fprintf('\n ----- Basic parameters to be "calibrated" -----');
fprintf('\n   Profit/value scaling      :     %4.0f', scaling);
fprintf('\n   Discount factor (beta)    :     %1.4f', beta);
fprintf('\n   Adoption cost drop (delta):     %1.4f', delta);
fprintf('\n   # potential entrants (Npe):     %1.0f', potent);
fprintf('\n -----------------------------------------------');
fprintf('\n');
fprintf('\n ------------ Estimated parameters -------------');
fprintf('\n   Operating fixed cost      :     %1.4f', phi);
fprintf('\n   Base SC for incumbents    :     %1.4f', kappa_inc);
fprintf('\n   Base SC for entrants      :     %1.4f', kappa_ent);
fprintf('\n -----------------------------------------------');
fprintf('\n');

% Profit scaling (Dan: "Treating data as parameters") & Value
Pi = Pi / scaling;  % Re-scaling profits, to avoid:
                    % (1) explosion of NaN's & (2) dominance of logit errors.
                    
% Terminal Values
V = zeros(T,4,12,12,15,12);     % Values (expected values, with expectation over shocks)
V(T,:,:,:,:,:) = Pi(T,:,:,:,:,:) / (1 - beta);  % V(T) pinned down by Pi(T)

% Initialize Solutions
EV = zeros(T,6,12,12,15,12);    % Expected values (expectation is over states)
Policy = zeros(T-1,6,12,12,15,12);% Choice probabilities

% Solve the model, just once
x0 = [phi kappa_inc kappa_ent];
LL = feval(@Likelihood, x0, 2);


%%%%%%%%%%%%%%%%%%%%%%%
% Post-solution tasks %
%%%%%%%%%%%%%%%%%%%%%%%
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task (1) Show the Values & Policies along the actual history (states). 

V_actual = zeros(18,4);
Policy_actual = zeros(17,6);

for t = 1:18
    
    No = State(t,1);
    Nb = State(t,2);
    Nn = State(t,3);
    Nroo = State(t,1);
    
    for type = 1:4
        V_actual(t,type) = V(t,type,No+1,Nb+1,Nn+1,Nroo+1);
    end
end

for t = 1:17
    
    No = State(t,1);
    Nb = State(t,2);
    Nn = State(t,3);
    Nroo = State(t,1);
    
    for type = 1:6
        Policy_actual(t,type) = Policy(t,type,No+1,Nb+1,Nn+1,Nroo+1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task (2) Simulate a history (based on the estimated model), and calc the "gap."

% Initialize
NS = 10000;                     % Number of simulations
npe = 4;                        % Number of potential entrants every year
State_simu = zeros(18,5,NS);    % Simulated history.
Adopt_simu = zeros(18,2,NS);    % Number of adoptions, by Inc & Ent.
seed = 645456;                  % Seed for random draws
rng(seed);                      % Set seed
U = rand(36+npe+11,17,NS);      % NS draws from U(0,1) to simulate each firm's discrete choice over 17 years
Uo = U(1:11,:,:);
Ub = U(12:22,:,:);
Un = U(23:36,:,:);
Upe= U(37:40,:,:);
Uroo=U(41:51,:,:);
Xo = zeros(11,17,NS);           % Number of exits among old-only firms
Eo = zeros(11,17,NS);           % Number of innovators among old-only firms
Xb = zeros(11,17,NS);           % Number of exits among both
Xn = zeros(14,17,NS);           % Number of exits among new-only
Epe= zeros(npe,17,NS);          % Number of innovators among potential entrants
Xroo=zeros(11,17,NS);           % Number of exits among incumbents' legacy divisions
xo = zeros(17,NS);
eo = zeros(17,NS);
xb = zeros(17,NS);
xn = zeros(17,NS);
epe= zeros(17,NS);
xroo=zeros(17,NS);

for ns = 1:NS                   % Loop 1: for each simulated history...
    fprintf('\n Simulating history #%5.0f', ns);
    State_simu(1,:,ns) = [State(1,:) State(1,1)];   % Year 1 is given by data
    State_simu(:,4,ns) = npe*ones(18,1);            % N_pe is fixed by assumption

    for t = 1:17                % Loop 2: for each year
        No = State_simu(t,1,ns);% Current state is given
        Nb = State_simu(t,2,ns);
        Nn = State_simu(t,3,ns);
        Npe = State_simu(t,4,ns);
        Nroo= State_simu(t,5,ns);
        PrXo = Policy(t,1,No+1,Nb+1,Nn+1,Nroo+1);
        PrEo = Policy(t,2,No+1,Nb+1,Nn+1,Nroo+1);
        PrXb = Policy(t,3,No+1,Nb+1,Nn+1,Nroo+1);
        PrXn = Policy(t,4,No+1,Nb+1,Nn+1,Nroo+1);
        PrEpe = Policy(t,5,No+1,Nb+1,Nn+1,Nroo+1);
        PrXroo = Policy(t,6,No+1,Nb+1,Nn+1,Nroo+1);
            
        if No > 0
            Xo(1:No,t,ns) = 1*(Uo(1:No,t,ns) < PrXo*ones(1:No,1));  % Old's exit decision
            xo(t,ns) = sum(Xo(1:No,t,ns),1);                        % Count of old's exits
            Eo(1:No,t,ns) = 1*(Uo(1:No,t,ns) > (1-PrEo)*ones(1:No,1));% Old's innovation decision
            eo(t,ns) = sum(Eo(1:No,t,ns),1);                        % Count of old's innovations
        else
            xo(t,ns) = 0;
            eo(t,ns) = 0;
        end
        
        if Nb > 0
            Xb(1:Nb,t,ns) = 1*(Ub(1:Nb,t,ns) < PrXb*ones(1:Nb,1));  % Both's exit decision
            xb(t,ns) = sum(Xb(1:Nb,t,ns),1);                        % Count of both's exits
        else
            xb(t,ns) = 0;
        end
        
        if Nn > 0
            Xn(1:Nn,t,ns) = 1*(Un(1:Nn,t,ns) < PrXn*ones(1:Nn,1));  % New's exit decision
            xn(t,ns) = sum(Xn(1:Nn,t,ns),1);                        % Count of new's exits
        else
            xn(t,ns) = 0;
        end
        
        if Npe > 0
            Epe(1:Npe,t,ns) = 1*(Upe(1:Npe,t,ns) < PrEpe*ones(1:Npe,1));% Entrants' innovation
            epe(t,ns) = sum(Epe(1:Npe,t,ns),1);                         % Count of entries
        else
            epe(t,ns) = 0;
        end
        
        if Nroo > 0
            Xroo(1:Nroo,t,ns) = 1*(Uroo(1:Nroo,t,ns) < PrXroo*ones(1:Nroo,1));  % Legacy's exit
            xroo(t,ns) = sum(Xroo(1:Nroo,t,ns),1);                              % Count of exits 
        else
            xroo(t,ns) = 0;
        end        
        
        State_simu(t+1,1,ns) = No - xo(t,ns) - eo(t,ns);    % No_prime
        State_simu(t+1,2,ns) = Nb - xb(t,ns) + eo(t,ns);    % Nb_prime
        State_simu(t+1,3,ns) = min([14, Nn - xn(t,ns) + epe(t,ns)]);    % Nn_prime (capped at 14)
        State_simu(t+1,5,ns) = Nroo - xroo(t,ns);           % Nroo_prime
        Adopt_simu(t+1,1:2,ns) = [eo(t,ns),epe(t,ns)];      % Counts of innovations
    end   
end

State_avg = sum(State_simu,3) / NS;     % Average over NS simulated histories
State_avg_rd = round(State_avg);        % Rounding firm counts to integers
Adopt_avg = sum(Adopt_simu,3) / NS;
Adopt_avg_rd = round(Adopt_avg);
   
% Calculate the "gap" index

Adopt_simu_rate = [Adopt_simu(:,1) / 11, Adopt_simu(:,2) / 14];
Adopt_simu_cdf = zeros(18,2);
for t = 2:18
    Adopt_simu_cdf(t,1) = Adopt_simu_cdf(t-1,1) + Adopt_simu_rate(t,1);
    Adopt_simu_cdf(t,2) = Adopt_simu_cdf(t-1,2) + Adopt_simu_rate(t,2);
end                
Adopt_simu_gap = Adopt_simu_cdf(:,2) - Adopt_simu_cdf(:,1);
Gap = sum(Adopt_simu_gap) * 100 / 17;   % In percentage, avg over 17 years.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task (3) Calculate the Social Welfare outcomes. ###

% Calculate CS (consumer surplus) in every year & industry state

fprintf('\n Making Consumer Surplus Matrix (T,No,Nb,Nn)=(18,12,12,15)...\n');
CS = zeros(18,12,12,15);  % (T x No x Nb x Nn)
for t = 1:18
    for n1 = 1:12
        for n2 = 1:12
            for n3 = 1:15
                
                % Adjustment: Matlab's index starts from 1
                no = n1 - 1;    
                nb = n2 - 1;
                nn = n3 - 1;
                
                % Prices are determined by the numbers of firms
                price_old = BigP(t,1,no+1,nb+1,nn+1);
                price_new = BigP(t,2,no+1,nb+1,nn+1);
                %price_old = MC(t,1) + (1 / (-alpha1 * (no + nb)));
                %price_new = MC(t,2) + (1 / (-alpha1 * (nn + nb)));
                
                % Deterministic parts of utilities in logit demand system
                delta_old = alpha0 + alpha1 * price_old + alpha2 * 0  ...
                            + alpha3 * X(t,1) + Xe(t,1) + Yd(t);
                delta_new = alpha0 + alpha1 * price_new + alpha2 * 1  ...
                            + alpha3 * X(t,2) + Xe(t,2) + Yd(t);
                
                % If a good doesn't exist, it's not in buyers' choice set
                if (no + nb) == 0; expDelta_old = 0;    % If Old goods don't exist
                else expDelta_old = exp(delta_old); end
                
                if (nn + nb) == 0; expDelta_new = 0;    % If New goods don't exist
                else expDelta_new = exp(delta_new); end
                
                % Consumer surplus = Market size * Individual CS
                CS(t,n1,n2,n3) = M(t) * (-1/alpha1) ...
                                * log(exp(0) + expDelta_old + expDelta_new);

%fprintf('\n Year %2.0f, State (n_old, n_both, n_new) = (%2.0f, %2.0f, %2.0f)',t,no,nb,nn);
%fprintf('\n     CS = %4.2f', CS(t,n1,n2,n3));
%fprintf('\n');

            end
        end
    end
end

% Calculate PS (producer surplus) in every year & industry state

fprintf('\n Making Producer Surplus Matrix (T,No,Nb,Nn)=(18,12,12,15)...\n');
PS = zeros(18,12,12,15);  % (T x No x Nb x Nn)
for t = 1:18
    for n1 = 1:12
        for n2 = 1:12
            for n3 = 1:15
                
                % Adjustment: Matlab's index starts from 1
                no = n1 - 1;    
                nb = n2 - 1;
                nn = n3 - 1;
                
                % Profits are determined by the numbers of firms
                pi_old  = Pi(t,1,n1,n2,n3);
                pi_both = Pi(t,2,n1,n2,n3);
                pi_new  = Pi(t,3,n1,n2,n3);
                
                % Producer surplus is just the aggregation of firm profits
                PS(t,n1,n2,n3) = no * pi_old + nb * pi_both + nn * pi_new;

%fprintf('\n Year %2.0f, State (n_old, n_both, n_new) = (%2.0f, %2.0f, %2.0f)',t,no,nb,nn);
%fprintf('\n     PS = %4.2f', PS(t,n1,n2,n3));
%fprintf('\n');

            end
        end
    end
end
PS = PS * 1000; % Re-scaling (b/c Pi is scaled down by 1000 when solving the model)

% Pick up CS & PS along the simulated history

fprintf('\n Picking CS & PS along the simulated history... \n');
CS_simu = zeros(18,1);
PS_simu = zeros(18,1);
for t = 1:18
    
    % Should count "Real Old-Only" as the producer of Old goods
    No = State_avg_rd(t,5);   
    Nb = State_avg_rd(t,2);
    Nn = State_avg_rd(t,3);
    
    CS_simu(t) = CS(t,No+1,Nb+1,Nn+1);
    PS_simu(t) = PS(t,No+1,Nb+1,Nn+1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the Dynamic Costs of Staying & Innovation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Counting the Dynamic Costs of Staying & Innovation... \n');

% Count the number of Exits along the simulated history
% (Adoptions are already counted in Counterfactual.m as Adopt_simu)
Adopt_simu2 = Adopt_avg_rd(2:18,1:2); % Watch out! - The timing of counting is different

for t = 1:17
        No = State_avg_rd(t,1);
        Nb = State_avg_rd(t,2);
        Nn = State_avg_rd(t,3);    
end
Nactive = sum(State_avg_rd(:,2:3),2) + State_avg_rd(:,5);   % Number of active firms

% Multiply # active firms & # adoptions by the cost parameters
FC_simu = phi * Nactive;     % Phi x # active firms
Delta = zeros(17,1);                % (17 x 1) vector of delta
for t = 1:17
    Delta(t) = delta^(t-1);
end
Kappa_inc = kappa_inc * Delta;       % Multiplier of Adoption Cost (inc)
Kappa_ent = kappa_ent * Delta;       % Multiplier of Adoption Cost (ent)
SC_simu = [Kappa_inc .* Adopt_simu2(:,1), Kappa_ent .* Adopt_simu2(:,2)];

% Rescale & collapse across types to summarize
SC_simu = SC_simu * 1000;
FC_simu  = FC_simu  * 1000;
SC_sum = sum(SC_simu,2);
FC_sum  = sum(FC_simu,2);
SC_sum = [SC_sum; 0];       % Just add an element for year 18 of Sunk Cost

% Social Welfare
SW = CS_simu + PS_simu + FC_sum - SC_sum;

% Get the discounted present value of the social welfare (as of 1981)
Beta = zeros(18,1); % Make a vector of the discount factor
for t = 1:18
    Beta(t,1) = beta ^ (t-1);
end
SWpv_1981thru98 = Beta' * SW;
SWpv_from1999   = beta ^ (19 - 1) * (SW / (1 - beta));
SWpv            = Beta .* SW;
SWpv            = [SWpv; SWpv_from1999];  
SWpv_sum        = SWpv_1981thru98 + SWpv_from1999;

% Also get the discounted present value of SW's components
CS_simu_pv = Beta' * CS_simu;     % PV of Consumer Surplus
PS_simu_pv = Beta' * PS_simu;     % PV of Producer Surplus
FC_simu_pv = Beta' * FC_sum;      % PV of Fixed Cost of Operation
SC_simu_pv = Beta' * SC_sum;      % PV of Sunk Cost of Innovation
SWpv_1981thru98_ByComponent = [CS_simu_pv, PS_simu_pv, FC_simu_pv, SC_simu_pv];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save & Finish
save Result_cf_NoCannibn.mat;
diary off;
warning on;