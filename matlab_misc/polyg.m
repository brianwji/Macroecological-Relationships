   function l = polyg(alpha , s);
        
     n = 1:1e7;
     l = sum( s .^ n ./ (n .^ alpha) );
        
    end