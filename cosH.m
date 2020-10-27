function C = cosH(Dataset, Scaling)
% Function to compute the hybrid-mahalanobis cosine between pairs of objects in an 
% n-by-m data matrix (n>>m).
%
% Input: 
% Dataset       a numerical data matrix of n rows (observations), and m variables 
%               or features with the condition n>>m.
%
% Scaling       "on" or "off" to scale the cosH scores between zero and 1;
%               default "on"
%
% Output:
%    C          pairwise cosine similarities, returned as a numeric 
%               row vector of length n(n?1)/2, corresponding to pairs of  
%               observations, where n is the number of observations in Dataset.
%               The distances are arranged in the order (1,2), (1,3), ..., (1,n),
%               (2,3), ..., (2,n), ..., (n-1,n), i.e., the upper triangle of 
%               the n-by-n similarity matrix.


% compute covariance
Cov = cov(Dataset);
n = size(Dataset,1);
m = size(Dataset,2);

if (n <= m)
    error(message('error: matrix dimensions not suitable: n<=m !'));
end

% compute eigenvalues
[coeff, eigValueDiag] = eig(Cov);
[eigValues, idx] = sort(diag(eigValueDiag), 'descend');
coeff = coeff(:,idx);
score = Dataset/coeff';

if any(eigValues<0)
    error(message('error: covariance matrix not positive semidefinite!'));
end

% compute pairwise cosine similarity
C = zeros(1,m*(m-1)/2);
d=0;
for i=1:n-1
   for j=i+1:n
   d=d+1;
   % inner product and norms
   inner = 0;
   normx1 = 0;
   normy1 = 0;
    
   for l=1:m        
  
       inner = inner + score(i,l).*score(j,l);  
       normx1 = normx1 + score(i,l).*score(i,l)./eigValues(l);
       normy1 = normy1 + score(j,l).*score(j,l)./eigValues(l);
 
   end
   
       C(d) = inner / (sqrt(normx1) * sqrt(normy1));    
      
   end
end

% normalization 
if (Scaling == "on")
C=C./max(abs(C));
end
end


