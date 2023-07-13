function  [H U]= jt_entropy(X,N)
%ENTROPY Compute the Shannon entropy of a set of variables.
%   ENTROPY(X,P) returns the (joint) entropy for the joint distribution 
%   corresponding to object matrix X and probability vector P.  Each row of
%   MxN matrix X is an N-dimensional object, and P is a length-M vector 
%   containing the corresponding probabilities.  Thus, the probability of 
%   object X(i,:) is P(i).  
%
if nargin <2
    N=100;
end

if length(unique(X))<500
    N = floor(sqrt(length(X))/2);
end

[nelements,centers]= hist(X,N);
nelements = nelements(find(nelements));
probSet = nelements./length(X);

%% Compute the entropy
H = -sum(probSet .* log2(probSet));   
%% compute the uniformity
U = sum(probSet .^2);  


