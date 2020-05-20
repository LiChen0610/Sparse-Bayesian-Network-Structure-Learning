function parents = SBN(X, lambda1, lambda2)
% This function implements the proposed algorithm.
% X is a n * p matrix (n: number of observation, p: number of nodes)
% lambda1, lambda2 correspond to symbol in the paper.
A = X;
[~,n] = size(A);
for i = 1:n
    t_mean = mean(A(:,i));
    t_std = std(A(:,i));
    A(:,i) =  (A(:,i) - t_mean)/t_std;
end
[~,prior_parents] = L1MB(A,n-1);
parents = prior_parents;
max_iter = 200;
iter = 0;
precision = 1e-6;
diff = 1;
while (iter <= max_iter) && (norm(diff,"fro") >= precision)
    parents_p = parents;
    iter = iter + 1;
    for j = 1:n
        dag = abs(parents) > 0;
        DG = mk_dag_dg(dag);
        [dist,~,~] = graphshortestpath(DG,j,'METHOD','BFS');
        dist = dist([1:j-1,j+1:end]);
        index = find(dist~=inf);
        path = zeros(n-1,1);
        path(index) = ones(length(index),1);
        w = SPDAS(A(:,[1:j-1,j+1:end]),A(:,j),parents([1:j-1,j+1:end],j),path,lambda1,lambda2);
        parents([1:j-1,j+1:end],j) = w;
    end
    diff = parents_p - parents;
end