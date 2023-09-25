function K = tvdlqr(A, B, Q, R, N)
    % A, B, Q, R are cell arrays with length N
    % A, B represent the state and input transition matrices, respectively
    % Q, R represent the state and input cost matrices, respectively

    % Initialize
    n = size(A{1},1); % dimension of the state
    P = cell(N,1);    % cell array to store cost-to-go
    P{N} = Q{N};      % final cost-to-go is final state cost

    % Compute the time-varying LQR controller gain K
    K = cell(N,1);
    for k = N:-1:2
        K{k-1} = (R{k-1} + B{k-1}' * P{k} * B{k-1}) \ (B{k-1}' * P{k} * A{k-1});
        P{k-1} = A{k-1}' * P{k} * A{k-1} + Q{k-1} + A{k-1}*P{k}* B{k-1}*K{k-1};
    end
end