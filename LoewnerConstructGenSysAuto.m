function [H] = LoewnerConstructGenSysAuto(SER,s)
% Get the inputs and outputs sizes
p =       SER.p;
m =      SER.m;

if p ~= m
    error("inputsN is not equal to outputsN")
end

% Number of ports is the same as inputsN if inputsN == outputsN
ports_N = p;
interp_N = length(s);   % Number of fine discretized samples

% Construct interpolant with generating system
% Construct Theta (and Theta_bar, replace ~ with <Theta_bar> for this option) 
[Theta,~] = LoewnerConstructGenSys(SER,s);

% Partition Block matrix function Theta into its blocks
Theta_11 = Theta(1:ports_N,1:ports_N,:);
Theta_12 = Theta(1:ports_N,ports_N+1:end,:);
Theta_21 = Theta(ports_N+1:end,1:ports_N,:);
Theta_22 = Theta(ports_N+1:end,ports_N+1:end,:);

% Allocate memory for the interpolant
H = zeros(ports_N,ports_N,interp_N);

G1 = eye(ports_N);
G2 = eye(ports_N);
    
Psi1 = pagemtimes(Theta_11,G1) - pagemtimes(Theta_12,G2);
Psi2 = pagemtimes(-Theta_21,G1) + pagemtimes(Theta_22,G2);

% Interpolant
H(:,:,:) = pagemrdivide(Psi1,Psi2);

