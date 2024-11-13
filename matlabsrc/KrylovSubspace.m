classdef KrylovSubspace
    properties
        A % Matrix A
        b % Vector b
        nb % Norm of b
        k % Krylov subspace dimension
        V % Krylov subspace basis
        H % Hessenberg or Tridiagonal matrix
        maxsize % Maximum size of the Krylov subspace
    end
    
    methods
        function obj = KrylovSubspace(A, b, maxk)
            % Initialize the Krylov subpace with the matrix A and the vector b and the maximum size of the Krylov subspace maxk
            obj.A = A;
            obj.b = b;
            obj.nb = norm(b,2);
            obj.maxsize = maxk;
            obj.k = 1;
            obj.V = zeros(length(b), maxk);
            obj.H = zeros(maxk+1, maxk+1);
            obj.V(:,1) = b/obj.nb;
        end
        
        function obj = arnoldi(obj)
            % Extend the Krylov subspace by one vector using the Arnoldi process
            if obj.k > obj.maxsize
                warning('Krylov subspace is full, you have to run either a restart or enlarge the subspace');
                return
            end
            w = obj.A*obj.V(:,obj.k);
            for j = 1:obj.k
                obj.H(j,obj.k) = obj.V(:,j)'*w;
                w = w - obj.H(j,obj.k)*obj.V(:,j);
            end
            obj.H(obj.k+1,obj.k) = norm(w,2);
            if obj.H(obj.k+1,obj.k) == 0
                warning('Krylov subspace has reached maximum dimension: lucky breakdown');
                return
            end
            obj.V(:,obj.k+1) = w/obj.H(obj.k+1,obj.k);
            obj.k = obj.k + 1;
        end

        function obj = reorth(obj)
            % Reorthogonalize the Krylov subspace basis using the modified Gram-Schmidt process
            for i = obj.k
                w = obj.V(:,i);
                for j = 1:i-1
                    obj.H(j,i) = obj.V(:,j)'*w;
                    w = w - obj.H(j,i)*obj.V(:,j);
                end
                obj.H(i+1,i) = norm(w,2);
                if obj.H(i+1,i) == 0
                    warning('Krylov subspace has reached maximum dimension: lucky breakdown');
                    return
                end
                obj.V(:,i+1) = w/obj.H(i+1,i);
            end
        end

        function obj = lanczos(obj)
            % Extend the Krylov subspace by one vector using the Lanczos process
            if obj.k > obj.maxsize
                warning('Krylov subspace is full, you have to run either a restart or enlarge the subspace');
                return
            end
            if obj.k == 1
                w = obj.A*obj.V(:,obj.k);
                alpha = obj.V(:,obj.k)'*w;
                w = w - obj.V(:,obj.k)*alpha;
                obj.H(obj.k,obj.k) = alpha;
                beta = norm(w,2);
                obj.V(:,obj.k+1) = w/beta;
                obj.H(obj.k+1,obj.k) = beta;
                obj.H(obj.k,obj.k+1) = beta;
            else
                w = obj.A*obj.V(:,obj.k) - obj.H(obj.k,obj.k-1)*obj.V(:,obj.k-1); 
                alpha = obj.V(:,obj.k)'*w;
                w = w - obj.V(:,obj.k)*alpha;
                obj.H(obj.k,obj.k) = alpha;
                beta = norm(w,2);
                obj.V(:,obj.k+1) = w/beta;
                obj.H(obj.k+1,obj.k) = beta;
                obj.H(obj.k,obj.k+1) = beta;
            end
            obj.k = obj.k + 1;
        end

        function obj = enlarge(obj,enlargefactor)
            % Enlarge the Krylov subspace by a factor of enlargefactor, i.e., new size = old size + round(old size*enlargefactor)
            
            %check if enlargefactor is a positive number
            if enlargefactor <= 0
                warning('Enlargement factor must be a positive number');
                return
            end

            obj.maxsize = obj.maxsize + round(obj.maxsize*enlargefactor);
            Vnew = zeros(length(obj.b), obj.maxsize);
            Hnew = zeros(obj.maxsize, obj.maxsize);
            Vnew(:,1:obj.k) = obj.V;
            Hnew(1:obj.k,1:obj.k) = obj.H;
            obj.V = Vnew;
            obj.H = Hnew;
        end

    end
end