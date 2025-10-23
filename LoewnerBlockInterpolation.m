function [SER] = LoewnerBlockInterpolation(H,s)
%
%   Standard block loewner 
%

% Assumes data on the form H \in \mathbb{C}^{p x m x N_s}

% Number of ports:
p = size(H,1);
m = size(H,2);

% Number of samples N_s
N_s = size(H,3);

% Rescale s
rescaling_factor = max(abs(s));

%rescalingFactor = 1;
s = s./rescaling_factor;
SER.rescaling_factor = rescaling_factor;

% Split data into mu and lambda
% First determine k and q
if mod(N_s,2) == 0
    % k = q+1, k+q <= N_s
    k = N_s/2 ;             %k is our split index, same notation as in literature
    q = N_s-k;

    % Alternate splitting
    indx_k = 1:2:N_s;
    indx_q = 2:2:N_s;

    % Left and Right splitting
    %indxk = 1:k;
    %indxq = k+1:N_s;

else
    % k = q+1, k+q <= N_s
    k =  (N_s+1)/2;            %k is our split index, same notation as in literature
    q = N_s-k;

    % Alternate splitting
    indx_k = 1:2:N_s;
    indx_q = 2:2:N_s;

    % Left and Right splitting
    %indxk = 1:k;
    %indxq = k+1:N_s;
end

mu = s(indx_q);
lambda = s(indx_k);

% Create Right data
Lambda = kron(diag(lambda),eye(m));
R = repmat(eye(m),[1, k]);
W = reshape(H(:,:,indx_k),[p, m*k]);

% Create Left data
M = kron(diag(mu),eye(p));
L = repmat(eye(p),[1,q]).';
V = reshape(pagetranspose(H(:,:,indx_q)),[m, p*q]).';

% Create Loewner matrix
temp1 = V*R - L*W;
temp2 = 1./(diag(M) - diag(Lambda).' );
Loewner_matrix = temp1.*temp2;

% Create shifted Loewner matrix for tangential data
temp3 = diag(M).*(V*R) - diag(Lambda).'.*(L*W);
Loewner_matrix_shifted = temp3 .* temp2;


% Assign E,A,B,C,D according to equation (11)
D = zeros(p,m);
L_right_inverse =Loewner_matrix\(V-L*D);   
A = Lambda + L_right_inverse*R;
B = L_right_inverse;
C = -(W - D*R);
E = eye(size(A));

% Save all matrices
model_order = size(A,1);
SER.p = p;
SER.m = m;
SER.E = E;
SER.A = A;
SER.B = B;
SER.C = C;
SER.D = D;
SER.Loewner_matrix_shifted = Loewner_matrix_shifted;    
SER.Loewner_matrix = Loewner_matrix;
SER.N_s = N_s;
SER.model_order = model_order;

% Add parameters to construct theta for Loewner with free parameter
SER.W = W;
SER.R = R;
SER.L = L;
SER.V = V;
SER.lambda = lambda;
SER.mu = mu;
SER.Lambda = Lambda;
SER.M = M;

end

