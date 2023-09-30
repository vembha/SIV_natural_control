#===========================
PLOTTER FILE FOR 
FIGURE S12B
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames, LsqFit, Distributions
using Statistics, DifferentialEquations, Random, Loess
Random.seed!(2023);
CairoMakie.activate!(type="svg")
size_cm = (5.0, 4.1);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the parameter files
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_1.xlsx","Sheet1";infer_eltypes=true))
df = df[!,2:13];

#=
Section 2: Define the required variables and obtain sampled parameter matrix
=#

n = 100_000; # Number of simulated individuals

# Define the matrices for storing results
mat_spVL  = Vector{Float64}[]; mat_tAUC = Vector{Float64}[];

#=
Section 3: Obtain the distributions for population parameters and sampling them
=#

# LogNormally distributed parameters
λ_dist_params = fit_mle(LogNormal, df.lambda_mode);
λ_dist = LogNormal(λ_dist_params.μ, λ_dist_params.σ);

dD_dist_params = fit_mle(LogNormal, df.dD_mode);
dD_dist = LogNormal(dD_dist_params.μ, dD_dist_params.σ);

r_dist_params = fit_mle(LogNormal, df.r_mode);
r_dist = LogNormal(r_dist_params.μ, r_dist_params.σ);

αE_dist_params = fit_mle(LogNormal, df.alphaE_mode);
αE_dist = LogNormal(αE_dist_params.μ, αE_dist_params.σ);

# Logit-Normally distributed parameters
fD_dist_params = fit_mle(LogitNormal, df.fD_mode);
fD_dist = LogitNormal(fD_dist_params.μ, fD_dist_params.σ);

# Normally distributed parameters
λE⁺_exp_dist_params = fit_mle(Normal, df.lambdaEstr_exp_mode);
λE⁺_exp_dist = Normal(λE⁺_exp_dist_params.μ, λE⁺_exp_dist_params.σ);

ω_exp_dist_params = fit_mle(Normal, df.omega_exp_mode);
ω_exp_dist = Normal(ω_exp_dist_params.μ, ω_exp_dist_params.σ);

# Fixed parameters
β_prime_exp = df.beta_prime_exp_mode[1]; μI = df.dI_mode[1];
θE = df.thetaE_mode[1]; μE = df.dE_mode[1]; T0_exp = df.T0_exp_mode[1];

# Parameter matrix given order: λ, β_prime_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
mat_param = hcat(
    rand(λ_dist, n), β_prime_exp.*ones(n, 1), rand(fD_dist, n), rand(λE⁺_exp_dist, n),
    μI.*ones(n, 1), rand(dD_dist, n), rand(r_dist, n), rand(αE_dist, n),
    θE.*ones(n, 1), μE.*ones(n, 1), rand(ω_exp_dist, n), T0_exp.*ones(n, 1)
);

#=
Section 3: Solving the model
=#

# Vectors to save the indices
idxs_SIC = Vector{Int64}[];
idxs_VIR = Vector{Int64}[];

# Vectors to save the virus dynamics for plot 5b-1
sol_V_SIC = Matrix{Float64}[];
sol_V_VIR = Matrix{Float64}[];

# Conditions required for ODE
time_vec = 0:0.1:365; tspan = (0, 365);
N = 1e6; C = 1e6;
sol_V = zeros(1,length(time_vec));
sol_ψ = zeros(1,length(time_vec));

# Best-fit model

for i = 1:1:n
    # Section 3a: Get the parameter and initial condition vectors
    curr_p    = Vector{Float64}(mat_param[i,:])
    init_vec  = [0.0,(10^(-2.76))/(curr_p[7]),0.0,0.0,0.0];
    curr_T0   = 10^(curr_p[12])
    curr_init = init_vec

    # Section 3b: Solve the model
    prob = ODEProblem(model_1!,curr_init,tspan,curr_p)
    soln = solve(prob,saveat = 0.1,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required solution vectors
    sol_V = log10.(2*curr_p[7].*max.(soln[2,:], vec((1e-8)*ones(1,length(time_vec)))))
    sol_ψ = soln[4,:].*soln[5,:]

    # Section 3d: Save the final estimates
    if sol_V[end] ≥ -3
        mat_spVL = vcat(mat_spVL, sol_V[end]) # Set-point viral load considered to be viral load at one year
        mat_tAUC = vcat(mat_tAUC, mid_point_int!(time_vec[1:281],sol_ψ[1:281])/28) # Time-averaged AUC of kE in the acute phase

        if sol_V[end] ≤ log10(400)
            if length(idxs_SIC) < 50 # 50 SIC dynamics are saved
                idxs_SIC  = vcat(idxs_SIC, i);
                if length(idxs_SIC) == 1
                    sol_V_SIC = vcat(sol_V_SIC, sol_V);
                else
                    sol_V_SIC = hcat(sol_V_SIC, sol_V);
                end
            end
        else
            if length(idxs_VIR) < 50 # 50 VIR dynamics are saved
                idxs_VIR  = vcat(idxs_VIR, i);
                if length(idxs_VIR) == 1
                    sol_V_VIR = vcat(sol_V_VIR, sol_V);
                else
                    sol_V_VIR = hcat(sol_V_VIR, sol_V);
                end
            end
        end
    end
end

mat_spVL  = convert(Vector{Float64}, mat_spVL);  mat_tAUC  = convert(Vector{Float64}, mat_tAUC);
idxs_SIC  = convert(Vector{Int64}, idxs_SIC);    idxs_VIR  = convert(Vector{Int64}, idxs_VIR);
sol_V_SIC = convert(Matrix{Float64}, sol_V_SIC); sol_V_VIR = convert(Matrix{Float64}, sol_V_VIR);

# Section 3e: Identifying SICs and VIRs

mat_spVL_SIC = Vector{Float64}[]; mat_spVL_VIR = Vector{Float64}[];
mat_tAUC_SIC = Vector{Float64}[]; mat_tAUC_VIR = Vector{Float64}[];
for i = 1:1:length(mat_spVL)
    if mat_spVL[i] <= log10(400)
        mat_spVL_SIC = vcat(mat_spVL_SIC, mat_spVL[i]);
        mat_tAUC_SIC = vcat(mat_tAUC_SIC, mat_tAUC[i]);
    else
        mat_spVL_VIR = vcat(mat_spVL_VIR, mat_spVL[i]);
        mat_tAUC_VIR = vcat(mat_tAUC_VIR, mat_tAUC[i]);
    end
end
mat_spVL_SIC = convert(Vector{Float64}, mat_spVL_SIC); mat_spVL_VIR = convert(Vector{Float64}, mat_spVL_VIR);
mat_tAUC_SIC = convert(Vector{Float64}, mat_tAUC_SIC); mat_tAUC_VIR = convert(Vector{Float64}, mat_tAUC_VIR);

#=
Section 4: Obtaining the fractions of SICs/VIRs in bins
=#

# Section 4a: Defining the required variables and the first-order Hill function model

b = 10; # Number of bins
thres_spVL = 200:50:600; # Different thresholds for spVL cutoff
thres_tAUC = zeros(length(thres_spVL), 1);
ζ_thres = 0.95;
@. sat(x, p) = p[1] + (1 - p[1])*x/(p[2] + x); # y = a + (1 - a)*x/(b + x)
p0 = [0.5, 0.5]; # Initial guesses for the parameters

# Section 4b: Estimating the fractions for each spVL threshold

for k = 1:1:length(thres_spVL) # For each spVL threshold
    curr_spVL = thres_spVL[k];
    SIC_fracs = zeros(b, 1); VIR_fracs = zeros(b, 1);
    for i = 1:1:b # For each bin
        n_SIC = 0; n_VIR = 0;
        for j = 1:1:length(mat_spVL) # For each element
            if mat_tAUC[j] > (i-1)/b && mat_tAUC[j] ≤ i/b
                if mat_spVL[j] ≤ log10(curr_spVL)
                    n_SIC = n_SIC + 1;
                else
                    n_VIR = n_VIR + 1;
                end
            end
        end
        SIC_fracs[i] = n_SIC/(n_SIC + n_VIR);
        VIR_fracs[i] = n_VIR/(n_SIC + n_VIR);
    end

    # Fitting the data to first-order Hill function
    fit_4f = curve_fit(sat, 0:0.1:0.9, vec(SIC_fracs), p0);
    a1 = coef(fit_4f)[1]; a2 = coef(fit_4f)[2];
    @. fit_model(x) = a1 + (1 - a1)*x/(a2 + x); # Final fit model
    thres_tAUC[k] = a2*((ζ_thres - a1)/(1 - ζ_thres));
end
thres_tAUC = vec(thres_tAUC);

#=
Section 5: Plot
=#

x_lim  = [180, 620];
y_lim  = [0, 1];
x_tick = [200, 400, 600];
y_tick = [0, 0.25, 0.5, 0.75, 1];

fig = Figure(resolution = size_pt, fontsize=12);

# Setting the panel
ax = Axis(fig[1, 1], spinewidth = 1, xticks = x_tick, yticks = y_tick,  xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

# Plotting commands
lines!(ax, thres_spVL, thres_tAUC, color = :steelblue3, linewidth = 2);
plot!(ax, thres_spVL, thres_tAUC, marker = :circle, strokecolor = :black, strokewidth = 1, color = :steelblue3, markersize = 6);

xlims!(ax, x_lim);
ylims!(ax, y_lim);

fig
save("fig_S12b.svg", fig, px_per_unit = 6)


