using  LinearAlgebra, Random, Distributions, Suppressor


######################################
#####Objective function  #############
######################################

function ϕ(t::Float64, gamma::Float64)
    return 1- exp(- t^2/gamma)
end


# esl empirical loss(l*)
function ESL(beta::Vector,dat::Array{Float64, 2}, gamma::Float64)
    n =size(dat)[1]
    X = dat[:,2:end]
    Y = dat[:,1]
    res = Y - X*beta
    Loss = sum( ϕ.(res, gamma))/n
    return Loss
end 

# adaptive lasso penalty
function adaLASSO(beta::Vector, lambda::Vector)
    return lambda'* abs.(beta)
end

# the objective function 
function ℓ(beta::Vector,dat::Array{Float64, 2}, gamma::Float64, lambda::Vector)
    return ESL(beta, dat, gamma) + adaLASSO(beta, lambda)
end


# true grad of l*
function ∇ESL( beta::Vector,dat::Array{Float64, 2}, gamma::Float64)
    X = dat[:, 2:end]
    Y = dat[:,1]
    n = size(dat)[1]
    res = X*beta - Y
    D = (res/gamma) .* exp.(-res.^2/gamma) .* X
    return  2 * mean(D, dims  = 1)[1,:]
end    
# we can also use ForwardDiff (autograd), which is equivalent



######################################
#####Proximal GD algorithm############
######################################

# soft-thresholding operator
function sfth(beta::Array{Float64,1}, lambda::Array{Float64,1})
    a = max.(abs.(beta) - lambda, 0 )
    return a.* sign.(beta)
end

# proximalGD algorithm of solving l1 penalized obj
function proxGD(x0::Array{Float64,1}, obj::Function, grd::Function, lrt::Function, sfth::Function, lambda::Array{Float64,1}, num_iters::Int = 100000)
    x = x0
    for i = 1:num_iters
        xx = x - lrt(i).* grd(x)
        x = sfth(xx, lrt(i)*lambda)
        if i%10000 == 0
            println(x, "\r")
        end
    end
    return x, obj(x)
end




##############################################
#### tuning parameter selection  #############
##############################################


# find pseudo outlier
function find_pseudo_out(beta::Vector,dat::Array{Float64, 2})
    X = dat[:, 2:end]
    Y = dat[:,1]
    n = size(dat)[1]
    res = X*beta - Y
    Sn = 1.4826 * median( abs.(res .- median(res)))
    ind = findall(x-> abs(x)≥ 2.5Sn, res)
    D_out = dat[ind, 1:end]
    D_good = dat[1:end .∉ [ind], 1:end]
    return D_out, D_good 
end

# use logdet for numerical stability
function logdet_V( beta::Vector,dat::Array{Float64, 2}, gamma::Float64)
    X = dat[:, 2:end]
    Y = dat[:,1]
    n = size(dat)[1]
    res = X*beta - Y
    #I^{-1}
    a = 2/gamma * mean(exp.(-res.^2/gamma).*(2*res.^2/gamma.-1))
    cov_x = X'*X/n
    inv_I = inv(cov_x)/a
    #Σ
    D = (2res/gamma) .* exp.(-res.^2/gamma) .* X
    st_D = D .- mean(D, dims = 1)
    Σ = st_D'*st_D/n
    return logdet(inv_I*Σ*inv_I)
end

#find range G
function ξ( beta::Vector,dat_good::Array{Float64, 2}, gamma::Float64, n)
    m = n - size(dat_good)[1]
    return 2m/n + 2/n*(n-m)*ESL(beta, dat_good, gamma)
end

function get_G(beta::Vector,dat_good::Array{Float64, 2}, gm::Vector, n)
    xxi = gm -> ξ(beta,dat_good, gm, n) 
    ξ_val = xxi.(gm)
    ind = findall(x -> x ≤ 1, ξ_val)
    G = gm[ind]
    return G, ξ_val[ind]
end

# find the optimal gamma
function search_gamma(beta::Vector, dat::Array{Float64, 2}, γ_seq::Vector) 
    n = size(dat)[1]
    _ ,dat_good = find_pseudo_out(beta, dat)  #find pseudo outliers
    xxi = gm -> ξ(beta,dat_good, gm, n)
    G, _ = get_G(beta, dat_good, γ_seq, n) # find G 
    dv = gm -> logdet_V(beta,dat, gm)
    gamma_opt = G[argmin(dv.(G))]  #find optimal γ
    return gamma_opt
end



# # making plots for the G and optimal gamma
# function gamma_plot(beta::Vector, dat::Array{Float64, 2}; gm::Vector= Array(.1:.1:20))
     # n = size(dat)[1]
    # _ ,dat_good = find_pseudo_out(beta, dat)  #find pseudo outliers
    
    # #p1 plot G and ξ val
    # xxi = gm -> ξ(beta,dat_good, gm, n)
    # p1 = plot(gm, xxi.(gm), ylim =(0,1.2) , lw = 3, xlabel = "γ", ylabel= "ξ(γ)" , title  = "Find G", label = "")
    # G, xi = get_G(beta, dat_good, Array(gm), n) # find G and corresponding ξ values
    # plot!(G, xi, lw=3, label = " G")
    
    # #p2 plot logdev and opt γ    
    # dv = gm -> logdet_V(beta,dat, gm ) 
    # p2 = plot(gm, dv.(gm), lw = 3, xlabel= "γ", ylabel= "logdet(V(γ))", title = "Find γ",label = "")
    # gamma_opt = search_gamma(beta, dat, gm)  #find optimal γ
    # plot!([gamma_opt], [minimum(dv.(gm))], mark = 5, label = "optimal γ")
    
    # # show the plot
    # plot(p1, p2, layout = (1,2))
# end





##########################################
### The complete procedure ###############
##########################################

function esl_lasso( beta0::Array{Float64,1}, dat::Array{Float64, 2}; γ_seq::Vector = Array(.1:.1:100),  opt_lrt::Function = i-> .1, opt_iters::Int = 50000, reps::Int = 2)
    ### x0: init_estimator (e.g. MM_estimator)
    ### dat: dataset s.t. Y = dat[:,1], X = dat[:, 2:end] 
    ### γ_seq: γ values to perform grid search (e.g. .1:.1:20)
    ### opt_alg: proxGD
    ### opt_lrt: smaller than Lipschitz
    ### opt_iters: niter for proxGD
    ### reps: number of reptitions of step1-3
    n =  size(dat)[1]
    X = dat[:,2:end]
    beta = beta0
    for k = 1:reps 
        # find the optimal gamma
        gamma_opt = search_gamma(beta, dat, γ_seq)
        println("γ = $gamma_opt")
        #adaLASSO parameter
        lambda = (log(n)/n)./ (abs.(beta) .+ 1e-10)  
        
        # update beta using proxGD
        obj = b-> ℓ(b, dat, gamma_opt, lambda)
        esl_grd = b-> ∇ESL(b, dat, gamma_opt)
        lip_cont = abs(eigmax(X'*X/n)*2/gamma_opt)
        lrt = i-> min(opt_lrt(i), lip_cont)
        beta, _ = proxGD(beta, obj, esl_grd, lrt, sfth, lambda, opt_iters)
    end
        
    return beta
end


function our_esl_lasso(beta0::Array{Float64,1}, dat::Array{Float64, 2}; opt_lrt::Function = i-> .1, opt_iters::Int = 50000, reps::Int = 2)
    n =  size(dat)[1]
    X = dat[:,2:end]
    beta = beta0
    
    for k  = 1:reps
        # find pseudo outlier 
        _, dat_good = find_pseudo_out(beta::Vector,dat::Array{Float64, 2})
        gamma =sum((dat_good[:,1] - dat_good[:,2:end]*beta).^2)/size(dat_good)[1]
        println("γ = $gamma")
        #adaLASSO parameter
        lambda = (log(n)/n)./ (abs.(beta) .+ 1e-10)  
        
        # update beta using proxGD
        obj = b-> ℓ(b, dat, gamma, lambda)
        esl_grd = b-> ∇ESL(b, dat, gamma)
        lip_cont = abs(eigmax(X'*X/n)*2/gamma)
        lrt = i-> min(opt_lrt(i), lip_cont)
        beta, _ = proxGD(beta, obj, esl_grd, lrt, sfth, lambda, opt_iters)
    end

    return beta
end

