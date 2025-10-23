function [indx_next, max_deviation, interpolants] = adaptiveSamplingTheta1(interp_data,freq_data,s_HQ,opts)
% Return estimated error also

% Set setting from options 
abs_or_rel = opts.abs_or_rel;
G_N = opts.G_N ;

% If predefined matrices G1 and G2 are used 
preDefinedGMat = 0;

% If predefined matrices are used
if isfield(opts,"GmatPredefined")
    preDefinedGMat = 1;
    G_N = size(opts.G1,3);
end

% Number of discrete samples 
interp_N = length(s_HQ);

% Perform interpolation
SER = LoewnerBlockInterpolation(interp_data,freq_data);

% Load additional parameters
p = SER.p;
m= SER.m;

% Construct Theta (and Theta_bar, replace ~ with <Theta_bar> for this option) 
[Theta,~] = LoewnerConstructGenSys(SER,s_HQ);

% Partition Block matrix function Theta into its blocks
Theta_11 = Theta(1:p,1:p,:);
Theta_12 = Theta(1:p,p+1:end,:);
Theta_21 = Theta(p+1:end,1:p,:);
Theta_22 = Theta(p+1:end,p+1:end,:);

if preDefinedGMat
    G1 = opts.G1;
    G2 = opts.G2;
else 
    G1 = 2*rand(p,m,G_N)-1;
    G2 = 2*rand(p,m,G_N)-1;
end

% Allocate memory for All of the interpolants
H_large = zeros(p,m,interp_N,G_N);

% Construct the interpolants
for ii = 1:G_N
    Psi1 = pagemtimes(Theta_11,G1(:,:,ii)) - pagemtimes(Theta_12,G2(:,:,ii));
    Psi2 = pagemtimes(-Theta_21,G1(:,:,ii)) + pagemtimes(Theta_22,G2(:,:,ii));

    H_large(:,:,:,ii) = pagemrdivide(Psi1,Psi2);
    %H_large(:,:,:,ii) = pagemtimes(Psi1,pagepinv(Psi2,1e-4));
end

% Calculate max deviation
if abs_or_rel  == "abs"

    max_deviation = 1/(p)*abs(range(H_large,4));
elseif abs_or_rel == "rel"

    % Technically an approximation here (not actually equation (25))
    max_deviation = 1/p*abs(range(H_large,4)./mean(H_large,4)); % Latest as of oct 8 2024

    % Exact equation (25) is as follows:
    % max_values = zeros(inputs_N,outputs_N,interp_N,G_N);
    % for ii = 1:G_N
    %     temp_max = max(abs((H_large - H_large(:,:,:,ii))./H_large),[],4);
    %     max_values(:,:,:,ii) = temp_max;
    % end

    % max_deviation2 = 1/inputs_N*max(max_values,[],[4]);

elseif abs_or_rel == "log"

    max_deviation = range(20*log10(abs(H_large)),4);
end


% Find maximum element entry
temp_max = max(max_deviation,[],[1 2],"linear");
[~,max_dev_indx] = max(temp_max);

% Select index corresponding to maximum element deviation as candidate
% frequency point
indx_next = max_dev_indx;

% Also return interpolants
interpolants = H_large;

end
