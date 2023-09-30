#==================================
FILE COMPRISING ALL THE ODE MODELS 
AS WELL AS THE FRAMEWORK
==================================#

#=
Please note that since death rate constants have the same notation 
as the differential (ex: dT == dT, dI == dI, and so on), μ has been 
used for the death rate constants to avoid confusion.
=#

function model_1!(du, u, p, t)

    T, I, D, E⁺, k⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, μE, ω_exp, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp); ω = 10^(ω_exp);

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I) - μE*E⁺
    du[5] = dk⁺ = ω*(1 - k⁺)
end

function model_1_sans_S!(du, u, p, t)

    T, I, D, E⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, μE, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp);

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I) - μE*E⁺
end

function model_2!(du, u, p, t)

    T, I, D, E⁺, k⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, αR, θX, μE, ω_exp, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp); ω = 10^(ω_exp); αX = αE + αR;

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I) - αX*E⁺*I/(θX + I) - μE*E⁺
    du[5] = dk⁺ = ω*(1 - k⁺)
end

function model_3!(du, u, p, t)

    T, I, D, E⁺, Q, k⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, ξ, ϕ, μE, κ, μQ, ω_exp, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp); ω = 10^(ω_exp); n = 1;

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I) - ξ*E⁺*((Q^n)/(ϕ^n + Q^n)) - μE*E⁺
    du[5] = dQ  = κ*I/(θE + I) - μQ*Q
    du[6] = dk⁺ = ω*(1 - k⁺)
end

function model_4!(du, u, p, t)

    T, I, D, E⁺, Q, k⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, ξ, ϕ, μE, κ, μQ, ω_exp, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp); ω = 10^(ω_exp); n = 4;

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I) - ξ*E⁺*((Q^n)/(ϕ^n + Q^n)) - μE*E⁺
    du[5] = dQ  = κ*I/(θE + I) - μQ*Q
    du[6] = dk⁺ = ω*(1 - k⁺)
end

function model_5!(du, u, p, t)

    T, I, D, E⁺, k⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, ρE, αE, θE, μE, ω_exp, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp); ω = 10^(ω_exp);

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + (ρE + αE*E⁺)*I/(θE + I) - μE*E⁺
    du[5] = dk⁺ = ω*(1 - k⁺)
end

function model_6!(du, u, p, t)

    T, I, D, E⁺, k⁺ = u; # Variables
    λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, μE, ω_exp, T0_exp = p; # Parameters
    β_prime = 10^(β_prime_exp); μT = (λ/(10^(T0_exp)));
    λE⁺ = 10^(λE⁺_exp); ω = 10^(ω_exp);

    du[1] = dT  = λ - β_prime*T*I - μT*T
    du[2] = dI  = (1 - fD)*β_prime*T*I - k⁺*E⁺*I - μI*I
    du[3] = dD  = fD*β_prime*T*I - μD*D
    du[4] = dE⁺ = λE⁺ + αE*E⁺*I/(θE + I + E⁺) - μE*E⁺
    du[5] = dk⁺ = ω*(1 - k⁺)
end

function S_func!(F)

    # Parameters for the ex vivo framework
    β = 1e-8; δ = 1.44; ρ = 0.36; p = 1440;
    c = 0.35; ϕ = p*ρ/c; f = 0.5; C = 1e6;
    T₀ = 1e6; V₀ = 10^2.86; τₘ = 5.64;

    # Framework
    γ₄ = δ;
    γ₈ = δ + F;

    α₄ = √((γ₄ - ρ)^2 + 4*f*β*T₀*ϕ);
    α₈ = √((γ₈ - ρ)^2 + 4*f*β*T₀*ϕ);

    V₄ = (0.5*V₀*exp(-0.5*τₘ*(γ₄ + ρ))/α₄)*((γ₄ - ρ + α₄)*(exp(-0.5*τₘ*α₄)) + (ρ - γ₄ + α₄)*(exp(0.5*τₘ*α₄)));
    V₈ = (0.5*V₀*exp(-0.5*τₘ*(γ₈ + ρ))/α₈)*((γ₈ - ρ + α₈)*(exp(-0.5*τₘ*α₈)) + (ρ - γ₈ + α₈)*(exp(0.5*τₘ*α₈)));

    return log10(V₄) - log10(V₈);
end

function mid_point_int!(x, y)
    I = 0
    for m = 2:length(x)
        I += 0.5*(y[m] + y[m-1])*(x[m] - x[m-1])
    end
    return I
end
