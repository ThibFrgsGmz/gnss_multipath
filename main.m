% [20/11/2017]
% Modelling the impact of Multipath on GNSS receivers
clear all;
clc;
close all;

% Simulation parameters 
min_corr = -1.5;  % minimum tau (chip) of the correlation function support
max_corr = 2;   % maximum tau (chip) of the correlation function support
tau_s    = 0.001; % sampling period (chip) of the correlation function

% Multipath envelope parameters
alpha    = 0.5;   % multipath relative amplitude
max_tau  = 1.5;   % maximum multipath relative delay
step_tau = 0.005;  % resolution of the envelope (chip)     
tau_vect = (0:step_tau:max_tau);
phi_vect = [0 pi]; % multipath phase offset relative to direct signal 

% Signal parameters
modulation = 'L1'; %
T_C        = 1; % chip duration (micro-seconds)

% Receiver parameters
c_s = 0.5; % receiver chip spacing (chip), in [0.05;1] range
B_W  = 1;  % receiver front-end bandwidth (MHz)
flag_double_delta = 0; % to be set to 1 if the double delta technique is used, 0 else

% Configuration parameters 
option.plot.ideal_corr_function       = 1; % 1 to plot correlation function of LOS
option.plot.mp_corr_function          = 0; % 1 to plot correlation function of LOS + Multipath
option.plot.filtered_corr_function    = 0; % 1 to plot correlation function of LOS + Multipath after front-end filtering
option.plot.envelope                  = 1; % 1 to plot multipath envelope
option.plot.position                  = 0; % 1 to plot the impact of the multipath in the position domain

% Generation of the correlation function of reference
if strcmp(modulation,'L1') == 1 % GPS L1 C/A
    
    tau = (min_corr:tau_s:max_corr); 
    K_0 = 1 - min(abs(tau),1); % autocorrelation of reference
    
    % plot ideal correlation function (no MP, no filtering)
    if option.plot.ideal_corr_function == 1  
        figure;
        plot(tau,K_0);
        title('Multipath-free correlation function');
        xlabel('code delay (chip)');
        ylabel('Correlation function');
        grid on;
        box on;
    end;
    
elseif strcmp(modulation,'E1') == 1 % GALILEO E1
    
    tau = (min_corr:tau_s:max_corr);
    
    %K_0 = ; % autocorrelation of reference, to be filled
    
    % plot ideal correlation function (no MP, no filtering)
    if option.plot.ideal_corr_function == 1  
        figure;
        plot(tau,K_0);
        title('Multipath-free correlation function');
        xlabel('code delay (chip)');
        ylabel('Correlation function');
        grid on;
        box on;
    end;
    
end;

for idx_phi = 1:length(phi_vect)
    
    phi_MP  = phi_vect(idx_phi);
    
    for idx_tau_MP = 1:length(tau_vect) % loop on the different possible multipath delays
        
        tau_MP = tau_vect(idx_tau_MP);
        
        K_mp  = K_0 + alpha * cos(phi_MP) * [zeros(1,round(tau_MP/tau_s)) K_0(1: end - round(tau_MP/tau_s))]; % correlation function with multipath      
       
        % plot correlation function with MP, no filtering (<=> infinite
        % front-end filtering)
        if option.plot.mp_corr_function == 1 && idx_tau_MP == floor(length(tau_vect)/2) && phi_MP == 0  
            figure;
            plot(tau,K_mp);
            title('Correlation function with multipath');
            xlabel('code delay (chip)');
            ylabel('Correlation function');
            grid on;
            box on;  
        end;
        
        % front-end filter modelling
        % generation of the filter
        fc12   = B_W/2;  % cutoff frequency (in Hz)
        Fe     = 1/(tau_s/T_C);
        nb_pts = length(K_mp);
        H      = zeros(1,nb_pts);
        H(1:round(fc12/Fe*nb_pts))                  = ones(1,round(fc12/Fe*nb_pts));
        H(nb_pts-round(fc12/Fe*nb_pts)+1:nb_pts)    = ones(1,round(fc12/Fe*nb_pts));
        F      = 0:Fe/nb_pts:Fe-Fe/nb_pts;
        % filtering by FFT (product in frequency domain)
        tf          = fft(K_mp);      % conversion in frequency domain
        tf          = tf.*H;          % application of filter
        K_filtered  = real(ifft(tf)); % conversion in time domain
        % end of front-end filtering

        % plot correlation function with MP, front-end filtering (<=> finite
        % front-end filtering)
        if option.plot.filtered_corr_function == 1 && idx_tau_MP == floor(length(tau_vect)/2) && phi_MP == 0
            figure;
            plot(tau,K_filtered);
            title('Correlation function with multipath filtered');
            xlabel('code delay (chip)');
            ylabel('Correlation function');
            grid on;
            box on; 
        end;

        if flag_double_delta == 0
            
            % create early and late correlation functions
            I_E = circshift(K_filtered, [0  -c_s/2/tau_s]);
            I_L = circshift(K_filtered, [0   c_s/2/tau_s]);
            
            % compute Early Minus Late discriminator
            Discriminator = I_E - I_L;
        
        elseif flag_double_delta == 1  % to be filled by students

            % create early and late and very early very late correlation functions

            % compute double delta discriminator

        end;
        
        % trick to limit the research
        st = find( tau(1:end)  < -c_s/2 ,1, 'last');
        en = find( tau(st:end) >  c_s/2 ,1, 'first');
        
        % find the zero crossing of the discriminator (the stability point)
        error = interp1(Discriminator(st:st+ en),tau(st:st+en) ,0,'linear');
        envelope(idx_phi,idx_tau_MP) = error ;
        
    end;
    
end;

if option.plot.envelope == 1
    % plot multipath envelope
    figure;
    plot(tau_vect,envelope(1,:),'b')
    hold on;
    plot(tau_vect,envelope(2,:),'b')
    grid on;
    xlabel('tau multipath (chip)');
    ylabel('DLL error (chip)');
end;




if option.plot.position  == 1
    
    MP_bias = max(max(abs(envelope)))*293; % conversion from chip to meters
    
    N = 1000;      % number of sample in the simulation
    sigma_UERE = 2; % sigma UERE in meters
    
    % Range domain to position domain
    H = [0.0225 0.9951 -0.0966 1;
        0.6750 -0.6900 -0.2612 1;
        0.0723 -0.6601 0.7477 1;
        -0.9398 0.2553 -0.2269 1;
        -0.5907 -0.7539 -0.2877 1;];
    
    % Measurement nominal errors generation
    Delta_Y = sigma_UERE*randn(5,N);
    
    % Add multipath bias on the chosen satellite
    for idx_t = 1:N
        Delta_X(:,idx_t) = inv(H'*H)*H'*Delta_Y(:,idx_t); %compute position error by LS
    end;
    
    figure;
    hold on;
    plot(Delta_X(1,:),Delta_X(2,:),'ob')
    plot(0,0,'or','markersize',6)
    legend('estimated position','true position')
    grid on;
    box on;
    xlabel('East (m)')
    ylabel('North (m)')
    
    figure;
    plot(sqrt(Delta_X(1,:).^2 + Delta_X(2,:).^2),'b')
    hold on;
    plot([1 N],mean(sqrt(Delta_X(1,:).^2 + Delta_X(2,:).^2))*[1 1],'r')
    grid on;
    legend('Position error','average Position error')
    box on;
    xlabel('time (s)')
    ylabel('Norm of horizontal position error (m)');
    
end;

