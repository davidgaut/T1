function [ys,MsD] = risky_steadystate(ys,M_,oo_,options_,nfun)
% David Gauthier - 2019 / Jan
% Risky steady state (De Groot 2013)
%-----------------------------------

%--------------------------------------------------------------------------------
% Moments must enter equations positively on the right hand side:
% y + x + M = 0;
% As in Dynare residuals for equation x = y + Mis computed as r = x - (y + M).
%--------------------------------------------------------------------------------

fun = str2func(nfun); if options_.order ~= 2; warning('Dynare approx order should be 2.'); end

% Provide Guess
improv = 1;
while improv > 1e-3
MsD_First = moment_compute(ys,M_,oo_,options_);
% ys_upd    = fun(MsD_First(options_.rss_id),ys,M_,options_.MPPrule);
ys_upd    = fun(MsD_First,ys,M_,options_.MPPrule);

improv    = norm(ys - ys_upd);
ys        = ys_upd;
end

% Solve
MsD = csolve(@tbn_fun,MsD_First(options_.rss_id),[],1e-14,5,ys,M_,oo_,options_,fun);
% MsD = csolve(@tbn_fun,MsD_First,[],1e-14,5,ys,M_,oo_,options_,fun);
ys  = fun(MsD,ys,M_,options_.MPPrule);
end


function [tbn,ys_upd,MsD] = tbn_fun(MsD,ys,M_,oo_,options_,fun)
    
    % Search Moment Convergence
    for i=1:size(MsD,2)
    MsD_try = MsD(:,i); 
    ys_upd  = fun(MsD_try,ys,M_,options_.MPPrule);

    try % Compute moments with last risky ss
        MsD_upd = moment_compute(ys_upd,M_,oo_,options_);
    catch
        tbn = inf(size(MsD,1),size(MsD,2)); 
        return
    end
    % Convergence of moments (solving eps_new = eps + Mom)
        tbn(:,i)  =  (MsD_try - MsD_upd(options_.rss_id));
%         tbn(:,i)  =  (MsD_try - MsD_upd);
    end
%     norm(tbn)

end

