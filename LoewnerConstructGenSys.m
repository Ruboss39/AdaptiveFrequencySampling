function [Theta,Theta_bar] = LoewnerConstructGenSys(SER,s)

% Load parameters
p =       SER.p;
m=      SER.m;
rescaling_factor = SER.rescaling_factor;
N_s = SER.N_s;
Loewner_matrix = SER.Loewner_matrix;
W = SER.W;
R = SER.R;
L = SER.L;
V = SER.V;
eye_pm = eye(p+m);
Lambda = SER.Lambda;
M = SER.M;

% Length of frequency parameter for constructing interpolant
s_N = length(s);

% -------------------------------------------------------------------------
% Correct the format of the frequency data, i.e., s and frequencyData, if necessary
% -------------------------------------------------------------------------

if size(s,1) > size(s,2)
    s = s.';
end

% Rescale s
s = s./rescaling_factor;

% Works for scalar as well as matrix case 
if mod(N_s,2) == 0
    % Factoring out Loewner matrix, calculation:
    s_diags = repmat(s,[N_s*p/2,1]);
    Lambda_diags = repmat(diag(Lambda),[1,s_N]);
    M_diags = repmat(diag(M),[1,s_N]);

    A1 = 1./(s_diags-Lambda_diags);
    A2 = 1./(s_diags-M_diags);

    X = Loewner_matrix\[L,V];
    Y = [-W;R]/Loewner_matrix;

    % dot multiplication of diag vector of A1 to X
    inverse_mat_1 = reshape(A1,[N_s*p/2,1,s_N]) .* X;
    inverse_mat_2 = reshape(A2,[1,N_s*p/2,s_N]) .* Y;

    Theta = eye_pm + pagemtimes([W;-R],inverse_mat_1);
    Theta_bar = eye_pm + pagemtimes(inverse_mat_2,[L,V]);
else
    % Factoring out Loewner matrix, calculation:
    s_diags1 = repmat(s,[(N_s+1)/2*p,1]);
    s_diags2 = repmat(s,[(N_s-1)/2*p,1]);

    Lambda_diags = repmat(diag(Lambda),[1,s_N]);
    M_diags = repmat(diag(M),[1,s_N]);

    A1 = 1./(s_diags1-Lambda_diags);
    A2 = 1./(s_diags2-M_diags);

    X = Loewner_matrix\[L,V];
    Y = [-W;R]/Loewner_matrix;

    % dot multiplication of diag vector of A1 to X
    inverse_mat_1 = reshape(A1,[(N_s+1)/2*p,1,s_N]) .* X;
    inverse_mat_2 = reshape(A2,[1,(N_s-1)/2*p,s_N]) .* Y;

    Theta = eye_pm + pagemtimes([W;-R],inverse_mat_1);
    Theta_bar = eye_pm + pagemtimes(inverse_mat_2,[L,V]);

end




    
