function y = zmatvec(A,D,x,transp)
%%ZMATVEC matrix vector product with Z = [O;I;I-D,A] or Z' = [O;I-D;I,A'] 

if ~exist("transp","var")
    transp = false;
elseif ~islogical(transp)
    transp = true;
end

n = size(A,1);
if transp
    % y = Z'*x
    y(1:n,1) = x(1:n,1)-D*x(1:n,1);
    y(n+1:2*n,1) = x(n+1:2*n,1)+A'*x(n+1:2*n,1);
else
    % y = Z*x
    y(1:n,1) = x(n+1:2*n,1);
    y(n+1:2*n,1) = x(1:n,1)-D*x(1:n,1)+A*x(n+1:2*n,1);
end

end


