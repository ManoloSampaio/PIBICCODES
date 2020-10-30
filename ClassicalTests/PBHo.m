function [eigenvectors,ranks] = PBHo(A,C)

[eigenvectors,eigenvalues] = eig(A); % Eigenvalues of A.
ranks=zeros(1,length(eigenvalues));  % Vector with the ranks

for i=1:length(eigenvalues)
    lambda = eigenvalues(i,i); % Set a eigenvalue for lambda.
    M = A-lambda*eye(size(A)); % Creates a matrix (A-I\lambda)
    ranked_matrix = [M;C]; 
    ranks(i)=rank(ranked_matrix);
end

end