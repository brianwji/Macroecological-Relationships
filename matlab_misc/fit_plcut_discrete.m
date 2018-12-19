function params = fit_plcut_discrete(data,xmin)

data = data(:);
N = numel(data);

% BPLCUTFIT fits a power-law cutoff model to the data
% Source: http://tuvalu.santafe.edu/~aaronc/powerlaws/bins/
%
%----------
% Options:
%----------
% 1. l = bplcutfit(h, boundaries, 'range', 1.5:0.01:3.5, 0.1:0.01:1)
%    The 'range' option can be specified to restrict search for
%    lambda and beta parameters. In above example, bplcutfit gives 
%    the best looking lambda and beta within the specified ranges. 
%    By default bplcutfit uses matlab's fminsearch function which 
%    in turn uses the Nelder-Mead simplex search algorithm. 
%    Refer following url:
%    (http://www.mathworks.com/help/techdoc/math/bsotu2d.html#bsgpq6p-11)
%
% 2. l = bplcutfit(h, boundaries, 'bmin', 100)
%    The 'bmin' option lets you fix a value such as 100 for bmin. 
%    Note that this value should be one of the values in the 
%    boundaries array. In the above example, 100 cannot be the 
%    last bin boundary. Also, it is advisable to give the fitting
%    procedure atleast two bins to work with. 
%    Note that if 'bmin' value is not one of the elements in 
%    boundaries, bplcutfit chooses the bin boundary which is closest 
%    to the specified value and less than that value. Also, the 
%    Default bmin value is the first bin boundary i.e. 
%    boundaries(1).
%
% Version 1.0 (2012)
% Copyright (C) 2012 Yogesh Virkar (University of Colorado, Boulder)
% Distributed under GNU GPL v3.0
% http://www.gnu.org/copyleft/gpl.html
% BPLCUTFIT comes with ABSOLUTELY NO WARRANTY


% Function to be minimized

    function l = polyg(alpha , s);
        
     n = 1:1e7;
     l = sum( s .^ n ./ (n .^ alpha) );
        
    end

function fval = plcutdiscrete_mle( vars )
    warning off all;
    
    alph = vars(1);
    lamb = vars(2);
    s = exp(-lamb);
    const_C = polyg( alph , s); %Normalization constant in PDF is 1/const_C
    
    % alpha should be greater than 1. 
    %if(vars(1)<1)
    %    fval = 10^10;
    %    return;
    %end
    
    fval =   -( -N * log(const_C) - alph * sum(log(data)) - lamb * sum(data));
    
    warning on all; 
end


% % Using Grid search 
%if ~isempty(rngal) && ~isempty(rnglam)
% rngal = [1:.1:2];
% rnglam = [.0001:.0001:.03];
%     logL = zeros(numel(rngal), numel(rnglam));
%     for i=1:numel(rngal)
%         for j=1:numel(rnglam)
%             logL(i,j) = plcutoff_mle([rngal(i) rnglam(j)]);
%             j
%         end
%         i
%     end
%     [C, I] = min(logL);
%     [~, I2] = min(C);
%     alpha_est = rngal(I(I2));
%     lambda_est = rnglam(I2);
    
% % Using fminsearch    
% else
    %Note that fminsearch could be sensitive to initial 
    %conditions. A reasonable initial condition could be 
    %startpt = [alpha, 0.0001], where alpha is the estimated alpha
    %value for power-law fit. 
    startpt = [1.9, 0.02];
    
    options.MaxFunEvals = 1e4;
    options.Display = 'iter';
    options.MaxIter = 1e4;
    vars_est = fminsearch(@plcutdiscrete_mle, startpt, options);
    alpha_est = vars_est(1);
    lambda_est = vars_est(2);
%end


L  = plcutdiscrete_mle([alpha_est,lambda_est]);

params.alpha_est = alpha_est;
params.lambda_est = lambda_est;
params.L = L;


end