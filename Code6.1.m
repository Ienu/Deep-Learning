visibleSize = 8*8;
hiddenSize = 25;
sparsityParam = 0.01;
lambda = 0.0001;
beta = 3;
%patches = sampleIMAGES;
b1 = zeros(hiddenSize, 1);
b2 = zeros(visibleSize, 1);
theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
%theta = initializeParameters(hiddenSize, visibleSize);
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
        t = 1;
    end
%     [t f g LSfunEvals] = WolfeLineSearch(x, t, d, f, g, gtd, 1e-4, 0.9, 4, 25, 1e-9, 0, 0, 1, @(p) sparseAutoencoderCost(p, ...
%         visibleSize, hiddenSize, ...
%         lambda, sparsityParam, ...
%         3, patches));
    x = x + t * d;
    [f g] = sparseAutoencoderCost(x, visibleSize, hiddenSize, lambda, sparsityParam, 3, patches);
end

W = reshape(x(1 : hiddenSize*visibleSize), hiddenSize, visibleSize);
display_network(W', 12);
