function y = phimatvec_mu(DL,A,D,Z,n,o,x,t,mu)
%%MATVEC Actual matrix-vector implementation for the exponential walk-based
%Laplacian. The quantity used in the routine are instantiated by the
%constructor buildphibasedlaplacian_mu.m. This function should be used through
%there.

if ~islogical(t)
    % PAY ATTENTION HERE: this is the converse of what it should be, is made
    % to be compatible with MATLABs expmv, and due to the fact that we wish
    % to solve the ODE with the transpose of this object.
    if strcmp('notransp',t)
        t = true;
    elseif strcmp('transp',t)
        t = false;
    elseif strcmp('real',t)
        y = 0;
        return 
    end
end

y = zeros(size(x));
if ~t % no transpose
    for i=1:size(x,2)
        tx = A*x(:,i);
        [v,stats] = phipm(1.0, Z, [o,[tx;A*tx - mu*D*x(:,i)]], 1e-12, false, 50);
        y(:,i) = DL.*x(:,i) - v(1:n,1);
    end
else % transpose 
    for i=1:size(x,2)
        [v,stats] = phipm(1.0, Z', [o,[o(1:n);x(:,i)]], 1e-12, false, 50);
        tv = A*v(1:n,1) + A*(A*v(1:n,1)) - mu*D*v(1:n,1);
        y(:,i) = DL.*x(:,i) - tv;
    end
end

fprintf('Number of substeps             : %d\n',stats(1));
fprintf('Number of rejected steps       : %d\n',stats(2));
fprintf('Number of Krylov steps         : %d\n',stats(3));
fprintf('Number of matrix exponentials  : %d\n',stats(4));

end