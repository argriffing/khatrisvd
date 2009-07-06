% add matlab7
% matlab -nojvm -nodisplay -nodesktop -nosplash -r svdtimer

function svdtimer()
n = 450;
p = 10000;
fprintf('create a %d by %d matrix:\n', n, p);
tic; X = rand(n, p); toc;
fprintf('take the svd of the matrix:\n');
tic; [U, S, V] = svd(X, 'econ'); toc;
quit;
