classdef MatrixFactorRegressor
    % MatrixFactorRegressor - Bilinear matrix regression for structured inputs
    %
    % Model: y ≈ tr(U' * A * V) for each matrix input A (size: m x k)
    %
    % Properties:
    %   U, V     - learned factor matrices
    %   rank     - low-rank dimension
    %   lambda   - regularization parameter
    %   max_iter - max number of iterations for training
    %
    % Methods:
    %   fit(A_all, y_all) - train the model
    %   predict(A_all)    - predict scalar outputs
    %
    
    properties
        U           % [m x r] factor
        V           % [k x r] factor
        rank        % rank r
        lambda      % regularization coefficient
        max_iter    % number of ALS iterations
    end
    
    methods
        function obj = MatrixFactorRegressor(r, lambda, max_iter)
            % Constructor
            if nargin < 3, max_iter = 50; end
            if nargin < 2, lambda = 1e-3; end
            obj.rank = r;
            obj.lambda = lambda;
            obj.max_iter = max_iter;
        end
        
        function obj = fit(obj, A_all, y_all)
            % Fit the model to training data
            n = length(y_all);
            [m, k] = size(A_all{1});
            r = obj.rank;
            lambda = obj.lambda;

            U = randn(m, r);
            V = randn(k, r);
            
            for iter = 1:obj.max_iter
                % Step 1: Fix V, update U
                ZU = zeros(size(U));
                for i = 1:n
                    ZU = ZU + -2/n*(y_all(i)-trace(U'*A_all{i}*V))*A_all{i}*V;
                end
                U = U + lambda * ZU;

                % Step 2: Fix U, update V
                ZV = zeros(size(V));
                for i = 1:n
                    ZV = ZV + -2/n*(y_all(i)-trace(U'*A_all{i}*V))*A_all{i}'*U;
                end
                V = V + lambda * ZV;
            end

            % Assign learned parameters
            obj.U = U;
            obj.V = V;
        end
        
        function y_pred = predict(obj, A_all)
            % Predict output for given input matrices
            n = length(A_all);
            y_pred = zeros(n, 1);
            for i = 1:n
                y_pred(i) = trace(obj.U' * A_all{i} * obj.V);
            end
        end
    end
end