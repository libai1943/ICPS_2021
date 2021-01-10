% Source codes for Cyber-Physical System based Cooperative Maneuver
% Planning and Replanning for Multiple Tractor-Trailer Vehicles in a
% Cluttered Environment. Bai Li 2021.
% Users must download AMPL together with IPOPT on their own.
clear all; close all; clc;

success_rate = [];
cpu_time = [];
num_iter = [];

for benchmark_case_id = 10 : 40
    clc
    global params_
    load(['benchmarks\', num2str(benchmark_case_id)]);
    %     params_.xmin = -20;
    %     params_.ymin = -20;
    %     params_.num_nodes_x = 40;
    %     params_.num_nodes_y = 40;
    %     params_.num_nodes_theta = 20;
    %     params_.resolution_x = 40 / params_.num_nodes_x;
    %     params_.resolution_y = 40 / params_.num_nodes_x;
    %     params_.resolution_theta = 2 * pi / params_.num_nodes_y;
    %     params_.M = [0 0 0];
    %     params_.R = 1.414;
    %     params_.Nobs = 15;
    %     params_.Nv = 2;
    %     params_.Obs = GenerateRandomObstacles();
    %     params_.costmap = CreateCostmap();
    %     params_.Nfe = 100;
    %     params_.config = GenerateVehicleConfig();
    params_.wheelbase = 1.5;
    params_.dt_in_hybrid_a_star_expansion = 2.0;
    params_.min_turning_radius = params_.wheelbase / tan(0.7);
    params_.penalty_multiplier_for_reversing = 3;
    params_.penalty_multiplier_for_drastic_direction_change = 5;
    params_.multiplier_H = 5;
    params_.phy_max = 0.7;
    params_.a_max = 0.5;
    params_.v_max = 3;
    params_.max_iter = 200;
    params_.num_iters_for_rs = 10;
    params_.terminal_xy_neiborhood = params_.dt_in_hybrid_a_star_expansion * 1;
    params_.terminal_theta_neiborhood = 1;
    WriteFilesForNLP();
    params_.terminal_xy_neiborhood = params_.dt_in_hybrid_a_star_expansion * 1;
    
    ig = cell(1, params_.Nv);
    for iv = 1 : params_.Nv
        [elem.x, elem.y, elem.theta] = PlanHybridAStarPath(iv);
        ig{1, iv} = elem;
        disp(['Warm start via FTMBHA algorithm for vehicle', num2str(iv)]);
    end
    
    [xc, yc] = WriteInitialGuess(ig);
    
    s_0 = 0;
    alpha = 2;
    beta = 0.5;
    max_iter = 500;
    maximum_allowable_s_ub = 3;
    
    % ACDO
    s_ub = s_0;
    WriteObstaclesForReducedNLP(xc, yc, s_ub);
    
    tic
    iter = 0;
    while (1)
        iter = iter + 1;
        !ampl r.run
        load sol_status.txt
        if (sol_status)
            [xc, yc] = LoadXY();
            if (IsSolValid(xc, yc))
                ready_flag = 1;
                break;
            end
            s_ub = s_ub + beta;
        else
            s_ub = s_ub - alpha;
        end
        s_ub = min(s_ub, maximum_allowable_s_ub);
        WriteObstaclesForReducedNLP(xc, yc, s_ub);
        if (iter > max_iter)
            ready_flag = 0;
            break;
        end
    end
    cpu_time = [cpu_time, toc];
    num_iter = [num_iter, iter];
    if (ready_flag)
        success_rate = [success_rate, 1];
        asd();
        drawnow;
    else
        success_rate = [success_rate, 0];
    end
end
disp(['Success rate = ', num2str(mean(success_rate) * 100), '%']);
disp(['Average CPU time = ', num2str(mean(cpu_time))]);
disp(['Average iteration numbers = ', num2str(mean(num_iter))]);