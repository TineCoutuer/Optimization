function grad = finite_diff(fun, x, idx)
    h = 1e-6;
    grad = zeros(1, length(x));
    for i = 1:length(x)
        x_fwd = x;
        x_bwd = x;
        x_fwd(i) = x_fwd(i) + h;
        x_bwd(i) = x_bwd(i) - h;

        f1 = fun(x_fwd);
        f2 = fun(x_bwd);

        if nargin == 3
            grad(i) = (f1(idx) - f2(idx)) / (2 * h);  % For specific constraint
        else
            grad(i) = (f1 - f2) / (2 * h);            % For objective
        end
    end
end
