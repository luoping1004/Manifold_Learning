% This script implements dgManifold which is an algorithm to predict
% disease-gene associations
clear;clc;

alpha = 0.1;
beta = alpha;
k = 100;

load('Adj.mat')
load('ModuleSim.mat')
load('Sg.mat')
load('Prior_wavg.mat')

[n,m] = size(adj);
distance = zeros(n,m);

Sd = Sd+eye(length(Sd));
Sg = Sg+eye(length(Sg));

A = adj;

%If the algorithm is trying to predict disease genes for a specific
%disease, add prior information for that disease. Add multiple rows of prior
%information when predicting associated genes for multiple diseases.

ind = 804; 

A(ind,:) = A(ind,:)+prior(ind,:);

Ar = diag(sum(A,2));
Ac = diag(sum(A));

Ld = diag(sum(Sd,2))-Sd;
Lg = diag(sum(Sg,2))-Sg;

L = [Ar+alpha*Ld -A; -A' Ac+beta*Lg];           

[V,D] = eig(L);

R = V(1:1611,2:k);
Q = V(1612:4977,2:k);

%Calculate the geodesic distance between all the genes and the target
%disease (Here the algorthm only calculates the distance for lung cancer)
i = ind;
for j=1:m
    if adj(i,j) == 1
        %Distance between the disease and known associated genes is set to
        %1e6
        distance(i,j) = 1e6;
    else
        distance(i,j) = sum((R(i,:)-Q(j,:)).^2);
    end
end
    
dname = ['distance_' str(ind) '.mat'];
save(dname, 'distance', '-mat');

