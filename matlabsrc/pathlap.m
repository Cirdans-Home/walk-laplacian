function Lpath = pathlap(G,type,parameter)
    %PATHLAP Builds the transformed k-Path Laplacian of the given type
    %   INPUT: G graph object
    %          type = 'exp' or 'power'
    %          parameter = scaling parameter
    %   OUTPUT: Lpath transformed path Laplacian
    
    Dist = G.distances();   % This may cost a lot!
    diam = max(Dist(:));
    Lpath = sparse(G.numnodes,G.numnodes);
    alpha = parameter;
    for i=1:diam
        [p,q] = find(Dist == i);
        e = ones(size(p));
        Gi = graph(sparse(p,q,e,G.numnodes,G.numnodes));
        Li = Gi.laplacian;
        if strcmpi(type,'exp')
            Lpath = Lpath + exp(-alpha*i)*Li;
        elseif strcmpi(type,'power')
            Lpath = Lpath + (1/i^alpha)*Li;
        else
            error('type has to be: exp or power');
        end
    end
    
    end
    
    