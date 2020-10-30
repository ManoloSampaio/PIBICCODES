%%%%
function [PBHtest,ranks] = PBHc(A,B)
dimension = length(A);
[eigenvectors,eigenvalues] = eig(A); % Eigenvalues of A.
ranks=zeros(1,length(eigenvalues));  % Vector with the ranks
PBHtest=1;
for i=1:length(eigenvalues)
    lambda = eigenvalues(i,i); % Set a eigenvalue for lambda.
    M = A-lambda*eye(size(A)); % Creates a matrix (A-I\lambda)
    ranked_matrix = [M B]; 
    ranks(i)=rank(ranked_matrix);
    if dimension~=rank(i)
        PBHtest=0; % If the rank condition is not respected, the test fails.
    end
end
    
end