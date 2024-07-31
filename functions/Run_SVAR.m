
%% Some Data Adjustments

y = data;

if ~isempty(exo)
    exog = y(:,exo);
    y    = y(:,endo);
end

if detrendIP,    y(:,2) = detrend(y(:,2));            end
if trend == 1,    exog = [(1:size(y,1))' exog];       end
if constant == 1, exog = [ones(size(y,1),1) exog];    end

[T,n] = size(y);
K     = n*p+size(exog,2);

dates_short = datetime(datevec(dates(p+1:end)),'Format','MM-yyyy');
T_short     = T - p;

maxtrys = 1e6;

%% Compute the closed-form solutions for the means
[~,~,posteriors,x,y_short] = BVAR(y,exog,p,prior_settings,1);
Beta_mean                  = reshape(posteriors.b_post,K,n);

%% Retrieve the S, Z, and Narrative matrices for identification

if ~exist('SR','var'),  SR = []; end
if ~exist('ZR','var'),  ZR = []; end
if ~exist('NSR','var'), NSR = []; end
if ~exist('EBR','var'), EBR = []; end

[S,Z,hirfs,Ns,Nc,~] = getRestrictions(SR,ZR,NSR,varNames,dates);

if  strcmp(StructuralIdentification,'Choleski')
    [BetaSigmaTries,Qs_per_BetaSigma] = deal(1);
end

if GLP == 1     % use Giannone et al (2015) to draw the hyperparameters
    Toutlier    = find(dates_short==datetime('07-2022','Format','MM-yyyy'));                                           % date of first outlier
    [b0,~,nu0,S0,Psi,pos,eta0] = get_prior_nc(y_short,n,p,K,prior_settings.hyperparameters,Toutlier,size(exog,2),sa_dummies,exog);         % priors & initial conditions
    [results,MIN,MAX]          = find_post_mode(prior_settings,y_short,x,p,T_short,n,K,b0,nu0,Psi,S0,eta0,pos,Toutlier,sa_dummies);       % VAR coefficients and residuals at the posterior model
    fprintf('Lambda_1 = %.2f, eta_1=%.2f, eta_2=%.2f, eta_2=%.2f, rho=%.2f\n',results.postmax.lambda1,results.postmax.eta(1),results.postmax.eta(2),results.postmax.eta(3),results.postmax.eta(4))
    P        = zeros(numDesiredDraws,length(results.postmode));                                                        % pre-allocation
    P(1,:)   = results.P;                                                                                              % use posterior mode as initial value
    HH       = results.HH;
    cholHH   = chol(HH);
    logMLold = -10e15;
    LOGML    = zeros(numDesiredDraws,1);
    while logMLold==-10e15          % get starting value for the Metropolis chain
        P(1,:)   = mvnrnd(results.P,HH,1,cholHH);
        logMLold = logMLVAR_formcmc(P(1,:)',prior_settings,y_short,x,p,T_short,n,K,b0,nu0,MIN,MAX,Psi,1,results.priorcoef,Toutlier,0,0,0,0,sa_dummies);
    end

    count   = 0;
    for i = 2:numDesiredDraws*2 % run the chain (we take some extra draws to ensure that we have enough)
        P(i,:)   = mvnrnd(P(i-1,:),HH,1,cholHH);
        logMLnew = logMLVAR_formcmc(P(i,:)',prior_settings,y_short,x,p,T_short,n,K,b0,nu0,MIN,MAX,Psi,1,results.priorcoef,Toutlier,0,0,0,0,sa_dummies);
        if logMLnew > logMLold                        % if likelihood increases, accept draw
            logMLold = logMLnew;
            count    = count + 1;
        else
            if rand(1) < exp(logMLnew-logMLold)       % accept also some "bad" draws
                logMLold = logMLnew;
                count    = count + 1;
            else                                      % draw is not accepted
                P(i,:)   = P(i-1,:);
                LOGML(i) = logMLold;
            end
        end
        if mod(i,500) == 0
            fprintf('Drawing hyperparameters, %d of %d\n',i,numDesiredDraws*2)
        end
    end
    fprintf('Hyperparameter estimation completed, acceptance ratio: %.1f%%\n',count/(numDesiredDraws*2)*100)
    if count/(numDesiredDraws*2)*100 < 15 || count/(numDesiredDraws*2)*100 > 35
        warning('acceptance ratio should be around 30%');
    end
end


z = 0;
%% Gibbs Sampler for Traditional Sign and Zero Restrictions
if tvp == 0
    S_post    = posteriors.S_post;
    v_post    = posteriors.v_post;
    b_post    = reshape(posteriors.b_post,K,n);
    CK        = chol(posteriors.NT,'lower')';
end

% Allocate Space for Saving Draws
Beta_save  = nan(K,n,numDesiredDraws);
Sigma_save = nan(n,n,numDesiredDraws);
A0_save    = nan(n,n,numDesiredDraws);
Q_save     = nan(n,n,numDesiredDraws);
w_save     = nan(1,numDesiredDraws);

numSavedDraws = 0;
numAttempts   = 0;

aux     = [zeros(n,n*p) ; eye(n*(p-1)) zeros(n*(p-1),n)];

while numSavedDraws < numDesiredDraws

    tic
    Beta_temp  = nan(K,n,BetaSigmaTries,Qs_per_BetaSigma);
    Sigma_temp = nan(n,n,BetaSigmaTries,Qs_per_BetaSigma);
    A0_temp    = nan(n,n,BetaSigmaTries,Qs_per_BetaSigma);
    Q_temp     = nan(n,n,BetaSigmaTries,Qs_per_BetaSigma);
    w_temp     = nan(1,BetaSigmaTries,Qs_per_BetaSigma);
    z          = z + 1;
    for worker = 1:BetaSigmaTries

        maxRoot = 2;
        trys    = 1;
        if GLP == 1  % draw coefficients conditional on estimated hyperparameters
            [B_draw,Sigma_draw,maxRoot,count_unstable] = draw_beta_sigma_GLP(P(z,:),y_short,x,p,T_short,n,K,b0,nu0,Psi,Toutlier,1,0,maxtrys,sa_dummies);
        else
            Sigma_draw = iwishrnd(S_post,v_post);
            COVARIANCE = chol(Sigma_draw,'lower');
            while maxRoot > .99 && trys <= maxtrys
                B_draw         = b_post + (CK\randn(K,n))*COVARIANCE';
                aux(1:n,1:n*p) = B_draw(size(exog,2)+1:end,:)';
                maxRoot        = max(abs(eig(aux)));
                trys           = trys + 1;
            end
        end

        if maxRoot <= .99 && ~isempty(B_draw)
            for Q_count = 1:Qs_per_BetaSigma
                if strcmp(StructuralIdentification,'Choleski')
                    Beta_temp(:,:,worker,Q_count)  = B_draw;
                    Sigma_temp(:,:,worker,Q_count) = Sigma_draw;
                    A0_temp(:,:,worker,Q_count)    = inv(chol(Sigma_draw));
                    Q_temp(:,:,worker,Q_count)     = eye(n);
                    w_temp(:,worker,Q_count)       = 1;
                else
                    % Check Sign and Zero Restrictions
                    [A0tilde_draw,Qdraw,uw] = SignAndZeroRestrictions(B_draw,Sigma_draw,p,exog,hirfs,S,Z,agnostic);
                    if isfinite(A0tilde_draw(1,1))
                        % Check Elasticity Bound Restrictions if Present
                        [checkElasticity] = ElasticityBoundRestrictions(EBR,varNames,shockNames,B_draw,A0tilde_draw,exog,n,p);
                        if checkElasticity
                            Beta_temp(:,:,worker,Q_count) = B_draw;
                            Sigma_temp(:,:,worker,Q_count) = Sigma_draw;
                            A0_temp(:,:,worker,Q_count) = A0tilde_draw;
                            Q_temp(:,:,worker,Q_count) = Qdraw;
                            w_temp(:,worker,Q_count) = uw;
                        end
                    end
                end
            end
        else
            %                 disp('unstable')
        end
    end
    Beta_clean  = reshape(Beta_temp,K,n,BetaSigmaTries*Qs_per_BetaSigma);
    Sigma_clean = reshape(Sigma_temp,n,n,BetaSigmaTries*Qs_per_BetaSigma);
    A0_clean    = reshape(A0_temp,n,n,BetaSigmaTries*Qs_per_BetaSigma);
    Q_clean     = reshape(Q_temp,n,n,BetaSigmaTries*Qs_per_BetaSigma);
    w_clean     = reshape(w_temp,1,BetaSigmaTries*Qs_per_BetaSigma);
    Beta_clean  = Beta_clean(:,:,isfinite(squeeze(A0_clean(1,1,:))));
    Sigma_clean = Sigma_clean(:,:,isfinite(squeeze(A0_clean(1,1,:))));
    A0_clean    = A0_clean(:,:,isfinite(squeeze(A0_clean(1,1,:))));
    Q_clean     = Q_clean(:,:,isfinite(squeeze(Q_clean(1,1,:))));
    w_clean     = w_clean(:,isfinite(w_clean));

    numAcceptedDraws = size(Beta_clean,3);

    if numAcceptedDraws > 0
        Beta_save(:,:,numSavedDraws+1:numSavedDraws+numAcceptedDraws)  = Beta_clean;
        Sigma_save(:,:,numSavedDraws+1:numSavedDraws+numAcceptedDraws) = Sigma_clean;
        A0_save(:,:,numSavedDraws+1:numSavedDraws+numAcceptedDraws)    = A0_clean;
        Q_save(:,:,numSavedDraws+1:numSavedDraws+numAcceptedDraws)     = Q_clean;
        w_save(:,numSavedDraws+1:numSavedDraws+numAcceptedDraws)       = w_clean;
        numSavedDraws                                                  = numSavedDraws + numAcceptedDraws;
    end

    numAttempts = numAttempts + (BetaSigmaTries*Qs_per_BetaSigma);

    disp(strcat('Attempted Draws:',num2str(numAttempts)))
    disp(strcat('Successful Draws:',num2str(numSavedDraws)))
    toc
end

%% Re-Weight Draws due to Zero Restrictions

if ~isempty(ZR)
    resample    = randsample(numSavedDraws,numSavedDraws,true,w_save);
    Beta_save   = Beta_save(:,:,resample);
    Sigma_save  = Sigma_save(:,:,resample);
    A0_save     = A0_save(:,:,resample);
    Q_save      = Q_save(:,:,resample);
end

%% Check Narrative Restrictions and Reweight

if ~isempty(NSR)
    disp('Checking Narrative Sign Restrictions')
    [Beta_narrative,Sigma_narrative,A0_narrative,weights_narrative] = CheckNarrativeAndReWeight(y,p,Beta_save,Sigma_save,A0_save,exog,Ns,Nc,nRepsWeights);
    numSavedNarrative                                               = size(Beta_narrative,3);
end

%% Compute IRFs and Structural Shocks for Traditional Sign Restrictions

if exist('numSavedDraws','var') == 0
    numSavedDraws = size(Beta_save,3);
end

if computeIRFs
    disp('Computing IRFs for traditional sign restrictions')
    hmax = 21;
    Draws_IRFs   = nan(n,n,hmax+1,numSavedDraws);
    Draws_Shocks = nan(T-p,n,numSavedDraws);
    for draw = 1:numSavedDraws
        IRFs                    = getIRFs(Beta_save(:,:,draw),A0_save(:,:,draw),exog,n,p,hmax);
        IRFs(cumulateWhich,:,:) = cumsum(IRFs(cumulateWhich,:,:), 3);
        Draws_IRFs(:,:,:,draw)  = IRFs;

        BETA = Beta_save(:,:,draw);
        B    = BETA(size(exog,2)+1:end,:)';
        c    = BETA(1:size(exog,2),:)';
        A0   = A0_save(:,:,draw);

        shocks                 = get_shocks(y,exog,A0,c,B,p);
        Draws_Shocks(:,:,draw) = shocks;
    end
else
    [Draws_IRFs,Draws_Shocks] = deal([]);
end

Draws_IRFs_const = Draws_IRFs;
Draws_Shocks_const = Draws_Shocks;


%% Compute IRFs and Structural Shocks that satisfy narrative restrictions
if computeIRFs

    if ~isempty(NSR)
        disp('Computing IRFs for narrative sign restrictions')
        numSavedNarrative      = size(Beta_narrative,3);
        Draws_IRFs_narrative   = nan(n,n,hmax+1,numSavedNarrative);
        Draws_Shocks_narrative = nan(T-p,n,numSavedNarrative);

        for draw = 1:numSavedNarrative
            IRFs                             = getIRFs(Beta_narrative(:,:,draw),A0_narrative(:,:,draw),exog,n,p,hmax);
            IRFs(cumulateWhich,:,:)          = cumsum(IRFs(cumulateWhich,:,:), 3);
            Draws_IRFs_narrative(:,:,:,draw) = IRFs;

            BETA = Beta_narrative(:,:,draw);
            B    = BETA(size(exog,2)+1:end,:)';
            c    = BETA(1:size(exog,2),:)';
            A0   = A0_narrative(:,:,draw);

            shocks                           = get_shocks(y,exog,A0,c,B,p);
            Draws_Shocks_narrative(:,:,draw) = shocks;
        end
    end
else
    [Draws_IRFs_narrative,Draws_Shocks_narrative] = deal([]);
end

Draws_IRFs_narrative_const   = Draws_IRFs_narrative;
Draws_Shocks_narrative_const = Draws_Shocks_narrative;


%% compute joint confidence sets for IRFs (Inoue/Kilian 2022)
% re-arrange draws
hmax = 21;
cl   = 0.68;             % confidence level for joint confidence set
IRF = []; IRF_narrative = [];

% one row for each model
for i = 1:n
    for j = 1:nshocks
        IRF           = [IRF squeeze(Draws_IRFs(i,j,1:hmax,:))'];
        IRF_narrative = [IRF_narrative squeeze(Draws_IRFs_narrative(i,j,1:hmax,:))'];
    end
end

% compute Bayes estimatior and joint confidence sets using absolute loss (Inoue/Kilian 2022)
[irfabs,jcs,I,S]                       = jointabsolute_c(IRF,cl);
[irfabs_narrative,jcs_narrative,In,Sn] = jointabsolute_c(IRF_narrative,cl);

% compute Bayes estimatior and joint confidence sets using quadratic loss (Inoue/Kilian 2022)
%     [irfabs,jcs,I]                         = jointquadratic_c(IRF,cl);
%     [irfabs_narrative,jcs_narrative,In,Sn] = jointquadratic_c(IRF_narrative,cl);


% re-arrange draws
Draws_JCS = nan(n,n,hmax,ceil(numSavedDraws*.68));      % holds joint confidence set
Draws_BE  = nan(n,n,hmax);                              % holds the bayes estimate

Draws_JCS_narrative = nan(n,n,hmax,ceil(numSavedNarrative*.68));
Draws_BE_narrative  = nan(n,n,hmax);

k = 1;
for i = 1:n
    for j = 1:nshocks
        Draws_BE(i,j,:)    = irfabs(1,(k-1)*(hmax)+1:k*(hmax));
        Draws_JCS(i,j,:,:) = jcs(:,(k-1)*(hmax)+1:k*(hmax))';

        Draws_BE_narrative(i,j,:)    = irfabs_narrative(1,(k-1)*(hmax)+1:k*(hmax));
        Draws_JCS_narrative(i,j,:,:) = jcs_narrative(:,(k-1)*(hmax)+1:k*(hmax))';
        k = k + 1;
    end
end

% compute simoultanous confidence bands (Montiel-Olea/Plaborg-Moller 2019)
Draws_IRFs_cum = permute(Draws_IRFs,[1,2,4,3]);
Draws_IRFs_cum = cumsum(Draws_IRFs_cum,3);

Draws_IRFs_narrative_cum = Draws_IRFs_narrative;
Draws_IRFs_narrative_cum = cumsum(Draws_IRFs_narrative_cum,3);
[Draws_Bands_narrative,Draws_Bands] = deal(nan(n,n,hmax,2));
for i = 1:n
    for j = 1:n
        Draws_Bands(i,j,:,:)           = SimInference.calibrated_Rbands(squeeze(Draws_IRFs(i,j,1:hmax,:))',cl)';                   % Calibrated Bayes
        Draws_Bands_narrative(i,j,:,:) = SimInference.calibrated_Rbands(squeeze(Draws_IRFs_narrative(i,j,1:hmax,:))',cl)';         % Calibrated Bayes
    end
end

