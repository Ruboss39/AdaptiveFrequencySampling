% example3.m             
%
% Adaptively sampling frequency samples from pre-calculated data.
% 
% This example script is part of the adaptive frequency sampling package
% (AdaptiveLoewnerTheta12.zip)
%
% Here, scalar data is interpolated.Additionally, all the interpolants are 
% plotted in each iteration
%
% Last revised: 17-01-2025 
% Created by:   Lucas Åkerstedt.
%

% Load example data 
% Scalar case
example_name = "4x1_Vivaldi_array";
ex_data_path = "./exampleData/4x1VivaldiArray.mat";
load(ex_data_path);
S_data = VivaldiArray_S;

% Testing scalar data
S_data = S_data(1,1,:);

% Sampling settings
adaptive_method = "Theta1";                  % "Theta1", "Theta2"
iterations = 22;                             % iterations to perform
cond_setting = 2;                            % How to calculate condition number [Theta2]
abs_or_rel = "rel";                          % rel, abs, log - minimize the relative, absolute or logarithmic error [Theta1]
G_N = 100;                                    % Number of random matrices to generate [Theta1] 
G_N_last = 600;                              % For estimating the final error (only on last round) [Theta1] 
double_side_sampling = 1;                    % Utilize passivity to sample H(s) and H(s)^*

% Plot all interpolants each iteration
plot_interpolants = 0;

% Interpolation settings
interp_data_complete = S_data;          % True data
freq_data_complete = 1i*reshape(frequency,1,[]);      % frequency parameter 1 x M_s
interp_N = length(freq_data_complete);  % maximum index
                 
s_HQ = freq_data_complete;              % Frequencies to create H(s_HQ) on (can be extended)

start_indx = [1 interp_N];              % Samples that we start with

frequency_samples_N = length(start_indx);

freq_indx = start_indx;

% Initial interpolation and frequency data
interp_data = interp_data_complete(:,:,freq_indx);
freq_data = freq_data_complete(freq_indx);

% Utilize passivity
if double_side_sampling
    freq_data = [freq_data, flip(conj(freq_data))];
    interp_data = cat(3,interp_data,flip(conj(interp_data),3));

end

% Previously selected points (their indices)
index_prev = freq_indx;

% Add options
opts.abs_or_rel = abs_or_rel;
opts.G_N = G_N;
opts.cond_setting = cond_setting;

% Iteration loop
for iter = 1:iterations 

    disp("Iteration: "+ string(iter));
    frequency_samples_N = frequency_samples_N+1;

    if adaptive_method == "Theta2"
        % Adaptive sampling
        [indx_next, max_deviation] = adaptiveSamplingTheta2(interp_data,freq_data,s_HQ,opts);

        % Add indx_next to index_prev
        index_prev = unique([index_prev, indx_next]);
    end
  

   if adaptive_method == "Theta1"
        % Adaptive sampling
        [indx_next, max_deviation, interpolants] = adaptiveSamplingTheta1(interp_data,freq_data,s_HQ,opts);
        %current_rel_error_F = squeeze(pagenorm(max_deviation,"fro"));
        %max_err = max(current_rel_error_F );
        %disp("Current max error: " + string(max_err));

        % Add indx_next to index_prev
        index_prev = unique([index_prev, indx_next]);
        
        % display interpolants
        if plot_interpolants
            p = size(interpolants,1);
            m = size(interpolants,2);
            interp_N = size(interpolants,3);
            G_N = size(interpolants,4);
            interps_dB = reshape(20*log10(abs(interpolants)),[p*m,interp_N,G_N]);
            sample_points_dB = reshape(20*log10(abs(interp_data)),p*m,[]);
            first_plot = 1;
            figInterp = figure(100);
            clf(figInterp);
            for ii = 1:G_N
                for jj = 1:p
                    plot(imag(s_HQ),interps_dB(jj,:,ii),'LineWidth',0.9);
                    if first_plot
                        hold on
                        first_plot = 0;
                    end
                    
                end
            end
            
            for jj = 1:p
                plot(imag(freq_data),sample_points_dB,'LineStyle','none','Marker','o','MarkerFaceColor','auto','MarkerEdgeColor','auto');
            end
    
            for jj = 1:p
                next_sample_dB = 20*log10(abs(reshape(interp_data_complete(:,:,indx_next),[],1)));
                plot(imag(s_HQ(indx_next)),next_sample_dB(jj),'LineStyle','none','Marker','*','MarkerFaceColor','auto','MarkerEdgeColor','auto');
            end
            
            xlim([min(imag(s_HQ)), max(imag(s_HQ))])
            grid on;
        end


   end

    % Update data
    freq_indx = [freq_indx,indx_next];

    % Sort indices
    [~, freq_indx_sorted] = sort(freq_indx);
    freq_indx = freq_indx(freq_indx_sorted);

    % Add new data
    interp_data = interp_data_complete(:,:,freq_indx);
    freq_data = freq_data_complete(freq_indx);

    % Utilize passivity
    if double_side_sampling
        freq_data = [freq_data, flip(conj(freq_data))];
        interp_data = cat(3,interp_data,flip(conj(interp_data),3));
    end

end

% Estimate error using Theta1 principle
opts.G_N = G_N_last;
[~, rel_error_est] = adaptiveSamplingTheta1(interp_data,freq_data,s_HQ,opts);


% Construct interpolant with sampled data
% Loewner interpolation
SER = LoewnerBlockInterpolation(interp_data,freq_data);

tic
% Construct interpolant using loewner generating system approach
[H2] = LoewnerConstructGenSysAuto(SER,s_HQ);
% Replace interpolation points with measured values 
H2(:,:,freq_indx) = interp_data_complete(:,:,freq_indx); 
InterpConstruct2 = toc

% Calculate absolute and relative error
abs_error2 = abs(H2-interp_data_complete);
rel_error2 = abs((H2-interp_data_complete)./interp_data_complete);

% Estimated error
rel_error_F = squeeze(pagenorm(rel_error2,"fro"));
rel_error_est = squeeze(pagenorm(rel_error_est,"fro"));

% Plot true data and interpolant
S_dB = 20*log10(abs(S_data));
H_dB = 20*log10(abs(H2));
S_dB_flat = reshape(S_dB,size(S_data,1)*size(S_data,2),[]);
H_dB_flat = reshape(H_dB,size(H2,1)*size(H2,2),[]);

figure();
for ii = 1:size(H2,1)*size(H2,2)
    plot(frequency,H_dB_flat,'LineWidth',0.9)
    if ii == 1
        hold on
    end
end
for ii = 1:size(H2,1)*size(H2,2)
    plot(frequency,S_dB_flat,'LineWidth',0.9,'Color','k','LineStyle','--')
end
xlim([min(frequency), max(frequency)])
ylim([-60 0])
xlabel('Frequency (GHz)')
ylabel('Magnitude (dB)')
grid on

figure();
semilogy(frequency,rel_error_F,'LineWidth',0.9);
hold on
semilogy(frequency,rel_error_est,'LineWidth',0.9,'Color','k','LineStyle','--');
xlim([min(frequency), max(frequency)])
xlabel('Frequency (GHz)')
ylabel('Magnitude (dB)')
xlim([min(frequency), max(frequency)])
ylim([1e-6, 1e0])
grid on
