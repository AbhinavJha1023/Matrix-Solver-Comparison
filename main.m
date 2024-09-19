timesGaussElimination = [];
timesGaussJordan = [];
timesGaussSeidel = [];
timesGaussJacobi = [];
timesCroutMethod = [];
timesDoolittleMethod = [];
fastestCounts = zeros(1, 6);
minTimesPerSize = inf(19, 6);
minMethodIndices = zeros(19, 1); 
for i = 2:20
    fastestCountCurrentSize = zeros(1, 6);
    for j = 1:50
        A = randi(100, i, i);
        B = randi(100, i, 1);
        tic;
        X = gaussElimination(A, B);
        timeGaussElimination = toc;
        timesGaussElimination = [timesGaussElimination, timeGaussElimination];
        tic;
        X = gaussJordan(A, B);
        timeGaussJordan = toc;
        timesGaussJordan = [timesGaussJordan, timeGaussJordan];
        tic;
        [X, iterations] = gaussSeidel(A, B, zeros(size(B)), 1e-3, 10000);
        timeGaussSeidel = toc;
        timesGaussSeidel = [timesGaussSeidel, timeGaussSeidel];
        tic;
        [X, iterations] = gaussJacobi(A, B, zeros(size(B)), 1e-3, 10000);
        timeGaussJacobi = toc;
        timesGaussJacobi = [timesGaussJacobi, timeGaussJacobi];
        tic;
        X = croutMethod(A, B);
        timeCroutMethod = toc;
        timesCroutMethod = [timesCroutMethod, timeCroutMethod];
        tic;
        X = doolittleMethod(A, B);
        timeDoolittleMethod = toc;
        timesDoolittleMethod = [timesDoolittleMethod, timeDoolittleMethod];
        times = [timeGaussElimination, timeGaussJordan, timeGaussSeidel, timeGaussJacobi, timeCroutMethod, timeDoolittleMethod];
        [minTime, fastestIndex] = min(times);
        fastestCounts(fastestIndex) = fastestCounts(fastestIndex) + 1;        
        if minTime < minTimesPerSize(i-1, fastestIndex) % i-1 for 1-based indexing
            minTimesPerSize(i-1, fastestIndex) = minTime;
            minMethodIndices(i-1) = fastestIndex;
        end
            fastestCountCurrentSize(fastestIndex) = fastestCountCurrentSize(fastestIndex) + 1;
    end
end
avgGaussElimination = mean(timesGaussElimination);
avgGaussJordan = mean(timesGaussJordan);
avgGaussSeidel = mean(timesGaussSeidel);
avgGaussJacobi = mean(timesGaussJacobi);
avgCroutMethod = mean(timesCroutMethod);
avgDoolittleMethod = mean(timesDoolittleMethod);
medianGaussElimination = median(timesGaussElimination);
medianGaussJordan = median(timesGaussJordan);
medianGaussSeidel = median(timesGaussSeidel);
medianGaussJacobi = median(timesGaussJacobi);
medianCroutMethod = median(timesCroutMethod);
medianDoolittleMethod = median(timesDoolittleMethod);
meanTimes = [avgGaussElimination, avgGaussJordan, avgGaussSeidel, avgGaussJacobi, avgCroutMethod, avgDoolittleMethod];
medianTimes = [medianGaussElimination, medianGaussJordan, medianGaussSeidel, medianGaussJacobi, medianCroutMethod, medianDoolittleMethod];
formattedMeanTimes = arrayfun(@(x) sprintf('%.7f', x), meanTimes, 'UniformOutput', false);
formattedMedianTimes = arrayfun(@(x) sprintf('%.7f', x), medianTimes, 'UniformOutput', false);
resultsMatrix = [formattedMeanTimes; formattedMedianTimes; num2cell(fastestCounts)];
methods = {'Gauss Elimination', 'Gauss-Jordan', 'Gauss-Seidel', 'Gauss-Jacobi', 'Crout Method', 'Doolittle Method'};
T = cell2table(resultsMatrix, 'VariableNames', methods, 'RowNames', {'Mean', 'Median', 'Fastest'});
ans1 = [];
for k = 1:length(minMethodIndices)
    minTime = minTimesPerSize(k, minMethodIndices(k));
    formattedMinTime = str2double(sprintf('%.7f', minTime));
    ans1 = [ans1, formattedMinTime];
    fprintf('For matrix size %d, the method with the least time was: %s (Time: %.7f)\n', k + 1, methods{minMethodIndices(k)}, minTimesPerSize(k, minMethodIndices(k)));
end
fprintf('\n');
disp(T);
[~, mostFrequentMethodIndex] = max(fastestCounts);
fprintf('The fastest method most frequently is: %s\n', methods{mostFrequentMethodIndex});
customIndices = 2:20;
figure;
bar(customIndices, ans1);
xlabel('Index');
ylabel('Value');
title('Bar Graph of Array Values vs. Custom Indices');
grid on;
function X = doolittleMethod(A, B)
    [n, ~] = size(A);
    L = eye(n); 
    U = zeros(n);
    for i = 1:n
        for j = i:n
            U(i,j) = A(i,j) - L(i,1:i-1) * U(1:i-1,j);
        end
        for j = i+1:n
            L(j,i) = (A(j,i) - L(j,1:i-1) * U(1:i-1,i)) / U(i,i);
        end
    end
    Y = zeros(n, size(B, 2));
    for i = 1:n
        Y(i, :) = (B(i, :) - L(i, 1:i-1) * Y(1:i-1, :));
    end
    X = zeros(n, size(B, 2));
    for i = n:-1:1
        X(i, :) = (Y(i, :) - U(i, i+1:n) * X(i+1:n, :)) / U(i, i);
    end
end
function X = croutMethod(A, B)
    n = length(B);
    L = zeros(n,n);
    U = eye(n);  
    for j = 1:n
        for i = j:n
            L(i,j) = A(i,j) - L(i,1:j-1)*U(1:j-1,j);
        end
        for k = j+1:n
            U(j,k) = (A(j,k) - L(j,1:j-1)*U(1:j-1,k)) / L(j,j);
        end
    end
    Y = zeros(n,1);
    for i = 1:n
        Y(i) = (B(i) - L(i,1:i-1)*Y(1:i-1)) / L(i,i);
    end
    X = zeros(n,1);
    for i = n:-1:1
        X(i) = (Y(i) - U(i,i+1:n)*X(i+1:n)) / U(i,i);
    end
end
function [X, iterations] = gaussJacobi(A, B, X, tolerance, max_iterations)
    n = length(B);
    iterations = 0;
    X_old = X;
    convergence = false;
    while ~convergence && iterations < max_iterations
        for i = 1:n
            sum = B(i);
            for j = 1:n
                if j ~= i
                    sum = sum - A(i, j) * X_old(j);
                end
            end
            X(i) = sum / A(i, i);
        end
        if norm(X - X_old, inf) < tolerance
            convergence = true;
        end
        X_old = X;
        iterations = iterations + 1;
    end
end
function [X, iterations] = gaussSeidel(A, B, X, tolerance, max_iterations)
    n = length(B);
    iterations = 0;
    convergence = false;
    while ~convergence && iterations < max_iterations
        X_old = X;
        for i = 1:n
            sum = B(i);
            for j = 1:n
                if j ~= i
                    sum = sum - A(i, j) * X(j);
                end
            end
            X(i) = sum / A(i, i);
        end
        if norm(X - X_old, inf) < tolerance
            convergence = true;
        end
        iterations = iterations + 1;
    end
end
function X = gaussJordan(A, B)
    [n, ~] = size(A);
    Aug = [A B];
    for i = 1:n
        [~, maxIndex] = max(abs(Aug(i:n, i)));
        maxIndex = maxIndex + i - 1;
        if maxIndex ~= i
            Aug([i, maxIndex], :) = Aug([maxIndex, i], :);
        end
        Aug(i, :) = Aug(i, :) / Aug(i, i);
        for j = 1:n
            if j ~= i
                factor = Aug(j, i);
                Aug(j, :) = Aug(j, :) - factor * Aug(i, :);
            end
        end
    end
    X = Aug(:, end);
end
function X = gaussElimination(A, B)
    [n, ~] = size(A);
    Aug = [A B];
    for i = 1:n
        [~, maxIndex] = max(abs(Aug(i:n, i)));
        maxIndex = maxIndex + i - 1;
        if maxIndex ~= i
            Aug([i, maxIndex], :) = Aug([maxIndex, i], :);
        end
        for j = i+1:n
            factor = Aug(j, i) / Aug(i, i);
            Aug(j, i:end) = Aug(j, i:end) - factor * Aug(i, i:end);
        end
    end
    X = zeros(n, size(B, 2));
    for i = n:-1:1
        X(i, :) = (Aug(i, end-size(B, 2)+1:end) - Aug(i, i+1:n) * X(i+1:n, :)) / Aug(i, i);
    end
end
