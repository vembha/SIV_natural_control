#===========================
PLOTTER FILE FOR 
FIGURES 3
REVISION VERSION #3
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
using LsqFit, HypothesisTests, Distributions
using Colors, ColorSchemes
CairoMakie.activate!(type="svg")
size_cm = (7.31, 5.42);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the data file
df_data = DataFrame(CSV.File("../../figure_data_files/v9_files/master_data_RNA_DNA_p27_kE0.csv"))

# Read the parameter file
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1.xlsx","Sheet1";infer_eltypes=true))

macaque_names = ["29915", "29925", "31041", "AV979",
"BA081", "BA209", "BB598", "BC094",
"BC179", "BC657", "BD536", "BD885",
"BG927", "BL669", "BO186", "BO413"]
# Progressors are placed in the top row, while controllers are placed in the remaining slots
ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]; # Order in which subplots must be placed

#=
Section 2: Define the required storage variables
=#

time_vec = 0:0.01:28; tspan = (0,28);
N = 1e6; C = 1e6;
sol_V = zeros(16,length(time_vec));
sol_I = zeros(16,length(time_vec));
sol_F = zeros(16,length(time_vec));
sol_S = zeros(16,length(time_vec));
sol_ψ = zeros(16,length(time_vec));

#=
Section 3: Solve the model for each macaque and store the values
=#

for i = 1:1:16
    # Section 3a: Get the parameter and initial condition vectors
    curr_p    = Vector{Float64}(df_param[i,2:13])
    init_vec  = [0.0,(10^(-2.76))/(curr_p[7]),0.0,0.0,0.0];
    curr_T0   = 10^(curr_p[12])
    curr_init = init_vec
    if i == 8 || i == 9 || i == 11 || i == 15
        curr_init[1] = curr_T0
    else
        curr_init[1] = curr_T0
        curr_init[2] = 10*curr_init[2]
    end

    # Section 3b: Solve the model
    prob = ODEProblem(model_1!,curr_init,tspan,curr_p)
    soln = solve(prob,saveat = 0.01,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required values and save them
    sol_V[i,:] .= log10.(2*curr_p[7].*max.(soln[2,:], vec((1e-8)*ones(1,length(time_vec)))))
    sol_I[i,:] .= log10.(soln[2,:] + soln[3,:])
    sol_F[i,:] .= (soln[4,:]./N).*C.*(soln[5,:])
    sol_ψ[i,:] .= (soln[4,:]./N)

    # Section 3d: Estimate S from the predictions for S
    sol_S[i,:] .= S_func!.(sol_F[i,:])
end

#=====================
FIGURE 3C:
S VS. TIME PREDICTED 
BY THE BEST-FIT MODEL 
FOR THE FIRST 28 DAYS
=====================#

#=
Section 4a: Plot
=#

x_lim  = [-1.5, 28];
y_lim  = [-0.1, 1.5];
x_tick = [0, 7, 14, 21, 28];
y_tick = [0, 0.5, 1.0, 1.5];

fig = Figure(resolution = size_pt, fontsize=12);
ax = Axis(fig[1,1], spinewidth = 1, xtickwidth = 1, ytickwidth = 1, xticks = x_tick, yticks = y_tick, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

# # Define the colors required for plotting
# IS_estimate = df_param.IS_mode; IS_estimate = vec(IS_estimate); IS_estimate = IS_estimate./maximum(IS_estimate);
# my_col_scheme = cgrad([:grey, :red], [0.0, 1.0]);

for i = 1:1:16
    if i == 3 || i == 4 || i == 7 || i == 12
        lines!(ax, time_vec, sol_S[i,:], color = :red, linewidth = 2);
    else
        lines!(ax, time_vec, sol_S[i,:], color = :grey, linewidth = 2);
    end
    # lines!(ax, time_vec, sol_S[i,:], color = my_col_scheme[IS_estimate[i]], linewidth = 1);
    xlims!(ax, x_lim);
    ylims!(ax, y_lim);
end

# cbar = Colorbar(fig[1,10]); cbar.ticks = ([0.0, 1.0], ["Controllers", "Progressors"])

fig
save("fig_3c.svg", fig, px_per_unit = 6)

#=====================
FIGURE 3A:
COMPARISON OF λE⁺ 
BETWEEN CONTROLLERS 
AND NON-CONTROLLERS
=====================#

#=
Section 4c: Read the required parameters and plot
=#

# Read the parameter files
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param_SIC = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1.xlsx","Controllers";infer_eltypes=true))
df_param_VIR = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1.xlsx","Progressors";infer_eltypes=true))

# Best-fit model
αE_SIC = df_param_SIC.alphaE_mode;
αE_VIR = df_param_VIR.alphaE_mode;
θE_SIC = df_param_SIC.thetaE_mode;
θE_VIR = df_param_VIR.thetaE_mode;
ω_SIC  = df_param_SIC.omega_exp_mode;
ω_VIR  = df_param_VIR.omega_exp_mode;

αE_SIC_mean = mean(αE_SIC); αE_VIR_mean = mean(αE_VIR);
θE_SIC_mean = mean(θE_SIC); θE_VIR_mean = mean(θE_VIR);
ω_SIC_mean = mean(ω_SIC); ω_VIR_mean = mean(ω_VIR);

# Estimating λEs

λ_SIC  = df_param_SIC.lambda_mode;
dI_SIC = df_param_SIC.dI_mode;
T0_SIC = broadcast(^, 10, df_param_SIC.T0_exp_mode);
dT_SIC = λ_SIC./T0_SIC;
β_prime_SIC = broadcast(^, 10, df_param_SIC.beta_prime_exp_mode);
dE_SIC = df_param_SIC.dE_mode;

λ_VIR  = df_param_VIR.lambda_mode;
dI_VIR = df_param_VIR.dI_mode;
T0_VIR = broadcast(^, 10, df_param_VIR.T0_exp_mode);
dT_VIR = λ_VIR./T0_VIR;
β_prime_VIR = broadcast(^, 10, df_param_VIR.beta_prime_exp_mode);
dE_VIR = df_param_VIR.dE_mode;

λE⁺_SIC = broadcast(^, 10, df_param_SIC.lambdaEstr_exp_mode);
λE⁺_VIR = broadcast(^, 10, df_param_VIR.lambdaEstr_exp_mode);

# Plot
size_cm = (5.42, 5.42);
size_pt = 38 .* size_cm;
fig = Figure(resolution = size_pt, fontsize=12);

# Horizontal jitter for barplot points
x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[1, 1], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), yticks = [0, 1, 2, 3], spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(λE⁺_SIC), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(λE⁺_VIR), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, λE⁺_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, λE⁺_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]); ylims!(ax, [0, 3]);
fig
save("fig_3a.svg", fig, px_per_unit = 6)

MannWhitneyUTest(λE⁺_SIC, λE⁺_VIR)

#=====================
FIGURE 3B:
COMPARISON OF αE 
BETWEEN CONTROLLERS 
AND NON-CONTROLLERS
=====================#

# Plot
size_cm = (5.42, 5.42);
size_pt = 38 .* size_cm;
fig = Figure(resolution = size_pt, fontsize=12);

# Horizontal jitter for barplot points
x1 = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 12));
x2 = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 4));
x1_mean = 0.75 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
x2_mean = 2.25 .+ vec(rand(Uniform(-0.25, 0.25), 1, 100));
ax = Axis(fig[1, 1], xticks = ([0.75, 2.25], ["Controllers", "Progressors"]), yticks = ([0.5, 0.6, 0.7, 0.8, 0.9]), spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
hidexdecorations!(ax, ticks = false, ticklabels = false);
hideydecorations!(ax, ticks = false, ticklabels = false);
hidespines!(ax, :t, :r);

scatter!(x1_mean, fill(median(αE_SIC), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x2_mean, fill(median(αE_VIR), 100), marker = :rect, color = :black, markersize = 6);
scatter!(x1, αE_SIC, strokecolor = :black, strokewidth = 1, color = :grey, markersize = 6);
scatter!(x2, αE_VIR, strokecolor = :black, strokewidth = 1, color = :red, markersize = 6);
xlims!(ax, [0, 3]); ylims!(ax, [0.48, 0.9]);
fig
save("fig_3b.svg", fig, px_per_unit = 6)

MannWhitneyUTest(αE_SIC, αE_VIR)

