function y = katzmatvec(DL,A,D,n,Li,alpha,x,t)
%KATZMATVEC Computation of the matrix-vector product with the Katz
%walk-based Laplacian.
%   Detailed explanation goes here

if ~islogical(t)
    % PAY ATTENTION HERE: this is the converse of what it should be, is made
    % to be compatible with MATLABs expmv, and due to the fact that we wish
    % to solve the ODE with the transpose of this object.
    if strcmp('notransp',t)
        t = true;
    elseif strcmp('transp',t)
        t = false;
    elseif strcmp('real',t)
        y = -1;
        return
    end
end

y = zeros(size(x));
flag = zeros(size(x,2),1);
if ~t
    matvec = @(x) x -alpha*(A*x) - alpha^2*(x - D*x);
    for i=1:size(x,2)
        [y(:,i),flag(i)] = pcg(matvec,x(:,i),1e-9,n,Li,Li');
        y(:,i) = (1-alpha)*(DL.*x(:,i) - y(:,i));
    end
else
    matvec = @(x) x -alpha*(A'*x) - alpha^2*(x - D*x);
    for i=1:size(x,2)
        [y(:,i),flag(i)] = pcg(matvec,x(:,i),1e-9,n,Li,Li');
        y(:,i) = (1-alpha)*(DL.*x(:,i) - y(:,i));
    end
end

if any(flag > 0) 
    warning('PCG Failure with flag %d',flag(flag > 0));
end