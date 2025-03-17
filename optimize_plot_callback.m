function stop = optimize_plot_callback(x, ~, state)
    global opt_path;
    
    if strcmp(state, 'iter')
        opt_path = [opt_path; x(:)']; % Store iteration values
    end
    
    stop = false; % Continue optimization
end
