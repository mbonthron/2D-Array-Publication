%% Clear Everything so there are no stragglers
function [data] = run_BB(data)
    %% Add the Paths to the Required Functions
    addpath('Barenblatt')
              
    data.A0 = zeros(2*data.N_modes,1);
    data.A0(1) = sqrt(12*(data.b/2)^2 - data.t^2)/sqrt(3);
    data.A0(2) = 1e-2;
    
    % Determine the coefficient matrix
    data = determine_coefficient_matrix(data);
        
    % Determine Ahat from A
    [A0hat] = determine_Ahat_from_A(data.A0,data);
    
    data.A0hat =  A0hat;
    
    %% Prepare for Time Integration
    data = determine_starting_vals(data);
    
    %% Run Initial Time Integration
    data.imposed = true;
    
    T_end = 1.25/data.alpha;
    tic
    [t1,Ahat1] = ode45(@(t,A) BB_single_arch_ODE(t,A,data),linspace(0,T_end,50000),data.A0hat);
    toc
    
    %% Recover A
    A1 = determine_A_from_Ahat(t1,Ahat1',data)';
    
    %% Determine the forces
    Q = determine_force(A1,data);
    
    first_positive = find(Q > 0,1);
    Q = Q(first_positive:end);
    t1 = t1(first_positive:end);
    A1 = A1(first_positive:end,:);
    
    last_positive = find(Q > 0,1,'last');
    
    A0 = A1(last_positive+1,:)';
    t0 = t1(last_positive+1);
    
    %% Switch to Non Imposed Displacement
    data.imposed = false;
    
    tic
    [t2,A2] = ode45(@(t,A) BB_single_arch_ODE(t,A,data),[0 0.1*T_end],A0);
    toc
    
    %% Combine the Results
    A = [A1(1:last_positive-1,:) ; A2];
    t = [t1(1:last_positive-1) ; t2+t1(last_positive)];
    
    %%
    figure(1); clf; hold on
    plot(t1,Q)
    scatter(t1(last_positive),Q(last_positive),100,'filled')
    
    figure(2); clf; hold on
    plot(t,A(:,1:3))
    plot(t1(last_positive)*[1 1],ylim,"k--","LineWidth",2)
    
    %% Save the Results
    file_string = "BB  " + data.file_string;
    close all
    save(file_string)
end