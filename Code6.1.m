visibleSize = 8*8;
hiddenSize = 25;
sparsityParam = 0.01;
lambda = 0.0001;
beta = 3;
patches = sampleIMAGES;
b1 = zeros(hiddenSize, 1);
b2 = zeros(visibleSize, 1);
c1 = 1e-4;
c2 = 0.9;
% theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
theta = initializeParameters(hiddenSize, visibleSize);
addpath minFunc2/

x = theta;
t = 1;
[f g] = sparseAutoencoderCost(x, visibleSize, hiddenSize, lambda, sparsityParam, beta, patches);

for i = 1 : 400
    i
    W = reshape(x(1 : 25 * 64), 25, 64);
    display_network(W', 12);
    if i == 1
        d = -g;
        old_dirs = zeros(length(g),0);
        old_stps = zeros(length(d),0);
        H = 1;
    else
        y = g - g_old;
        s = t * d;
        ys = y' * s;
        
        numCorrections = size(old_dirs, 2);
        if numCorrections < 100
            old_dirs(:, numCorrections + 1) = s;
            old_stps(:, numCorrections + 1) = y;
        else
            old_dirs = [old_dirs(:, 2 : 100) s];
            old_stps = [old_stps(:, 2 : 100) y];
        end
        H = ys / (y' * y);
        [p k] = size(old_dirs);
        rou = 1 ./ sum(old_stps .* old_dirs);
        
        q = zeros(p, k + 1);
        r = zeros(p, k + 1);
        alpha = zeros(k, 1);
        beta = zeros(k, 1);
        
        q(:, k + 1) = -g;
        for ii = k : -1 : 1
            alpha(ii) = rou(ii) * old_dirs(:, ii)' * q(:, ii + 1);
            q(:, ii) = q(:, ii + 1) - alpha(ii) * old_stps(:, ii);
        end
        r(:, 1) = H * q(:, 1);
        
        for ii = 1 : k
            beta(ii) = rou(ii) * old_stps(: , ii)' * r(:, ii);
            r(:, ii + 1) = r(:, ii) + old_dirs(:, ii) * (alpha(ii) - beta(ii));
        end
        d = r(:, k + 1);
    end
    g_old = g;
    gtd = g' * d;
    if i == 1
        t = min(1, 1 / sum(abs(g)));
    else
        %% WofeLineSearch
        [f_new g_new] = sparseAutoencoderCost(x + t * d, visibleSize, hiddenSize, lambda, sparsityParam, 3, patches);
        gtd_new  = g_new' * d;
        LSiter = 0;
        t_prev = 0;
        f_prev = f;
        g_prev = g;
        gtd_prev = gtd;
        done = 0;
        while LSiter < 25
            if f_new > f + c1 * t * gtd || (LSiter > 1 && f_new >= f_prev)
                bracket = [t_prev t];
                bracketF = [f_prev f_new];
                bracketG = [g_prev g_new];
                break;
            elseif abs(gtd_new) <= -c2 * gtd
                bracket = t;
                bracketF = f_new;
                bracketG = g_new;
                done = 1;
                break;
            elseif gtd_new >= 0
                bracket = [t_prev t];
                bracketF = [f_prev f_new];
                bracketG = [g_prev g_new];
                break;
            end
            temp = t_prev;
            t_prev = t;
            minStep = t + 0.01 * (t - temp);
            maxStep = t * 10;
            
            A = [temp^3 temp^2 temp 1;
                t^3 t^2 t 1;
                3 * temp^2 2 * temp 1 0;
                3 * t^2 2 * t 1 0];
            b = [f_prev; f_new; gtd_prev; gtd_new];
            
            params = A \ b;
            
            dParams = zeros(3, 1);
            
            dParams(1) = params(1) * 3;
            dParams(2) = params(2) * 2;
            dParams(3) = params(3);
            
            if any(isinf(dParams))
                cp = [minStep; maxStep; temp; t].';
            else
                cp = [minStep; maxStep; temp; t; roots(dParams)].';
            end
            
            fmin = inf;
            t = (minStep + maxStep) / 2;
            for xCP = cp
                if imag(xCP) == 0 && xCP >= minStep && xCP <= maxStep
                    fCP = polyval(params, xCP);
                    if imag(fCP)==0 && fCP < fmin
                        t = real(xCP);
                        fmin = real(fCP);
                    end
                end
            end
            
            f_prev = f_new;
            g_prev = g_new;
            gtd_prev = gtd_new;
            [f_new g_new] = sparseAutoencoderCost(x + t * d, visibleSize, hiddenSize, lambda, sparsityParam, 3, patches);
            gtd_new = g_new' * d;
            LSiter = LSiter + 1;
        end
        
        if LSiter == 25
            bracket = [0 t];
            bracketFval = [f f_new];
            bracketGval = [g g_new];
        end
        
        insufProgress = 0;
        while ~done && LSiter < 25
            [f_LO LOpos] = min(bracketF);
            HIpos = 3 - LOpos;
            
            xmin = min(bracket(1), bracket(2));
            xmax = max(bracket(1), bracket(2));
            
            A = [bracket(1)^3 bracket(1)^2 bracket(1) 1;
                bracket(2)^3 bracket(2)^2 bracket(2) 1;
                3 * bracket(1)^2 2 * bracket(1) 1 0;
                3 * bracket(2)^2 2 * bracket(1) 1 0];
            
            b = [bracketF(1); bracketF(2); bracketG(:, 1)' * d; bracketG(:, 2)' * d];
            
            params = A \ b;
            
            dParams = zeros(3, 1);
            
            dParams(1) = params(1) * 3;
            dParams(2) = params(2) * 2;
            dParams(3) = params(3);
            
            if any(isinf(dParams))
                cp = [xmin; xmax; bracket(1); bracket(2)].';
            else
                cp = [xmin; xmax; bracket(1); bracket(2); roots(dParams)].';
            end
            
            fmin = inf;
            t = (xmin + xmax) / 2;
            for xCP = cp
                if imag(xCP) == 0 && xCP >= xmin && xCP <= xmax
                    fCP = polyval(params, xCP);
                    if imag(fCP)==0 && fCP < fmin
                        t = real(xCP);
                        fmin = real(fCP);
                    end
                end
            end

            if min(max(bracket) - t, t - min(bracket)) / (max(bracket) - min(bracket)) < 0.1
                if insufProgress || t >= max(bracket) || t <= min(bracket)
                    if abs(t - max(bracket)) < abs(t - min(bracket))
                        t = max(bracket) - 0.1 * (max(bracket) - min(bracket));
                    else
                        t = min(bracket) + 0.1 * (max(bracket) - min(bracket));
                    end
                    insufProgress = 0;
                else
                    insufProgress = 1;
                end
            else
                insufProgress = 0;
            end
            
            [f_new g_new] = sparseAutoencoderCost(x + t * d, visibleSize, hiddenSize, lambda, sparsityParam, 3, patches);
            gtd_new = g_new' * d;
            LSiter = LSiter + 1;
            
            if f_new > f + c1 * t * gtd || f_new >= f_LO
                bracket(HIpos) = t;
                bracketF(HIpos) = f_new;
                bracketG(:, HIpos) = g_new;
            else
                if abs(gtd_new) <= -c2 * gtd
                    done = 1;
                elseif gtd_new * (bracket(HIpos) - bracket(LOpos)) >= 0
                    bracket(HIpos) = bracket(LOpos);
                    bracketF(HIpos) = bracketF(LOpos);
                    bracketG(:, HIpos) = bracketG(:, LOpos);                    
                end
                bracket(LOpos) = t;
                bracketF(LOpos) = f_new;
                bracketG(:, LOpos) = g_new;
            end
            
            if ~done && abs((bracket(1) - bracket(2)) * gtd_new) < 1e-9
                break;
            end
        end
        
        [f_LO LOpos] = min(bracketF);
        t = bracket(LOpos);
    end
    x = x + t * d;
    [f g] = sparseAutoencoderCost(x, visibleSize, hiddenSize, lambda, sparsityParam, 3, patches);
end

W = reshape(x(1 : hiddenSize*visibleSize), hiddenSize, visibleSize);
display_network(W', 12);
