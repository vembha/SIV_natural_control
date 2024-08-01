#===========================
PLOTTER FILE FOR 
FIGURE S15A
REVISION VERSION #3
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
CairoMakie.activate!(type="svg")
size_cm = (18, 18);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the data file
df_data = DataFrame(CSV.File("../../figure_data_files/v9_files/master_data_RNA_DNA_p27_kE0.csv"))

# Read the parameter files
# Parameter order (Total 12): λ, β_hat_exp, fD, λE⁺_exp, dI, dD, r, αE, θE, dE, ω_exp, T0_exp
df_param       = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1.xlsx","Sheet1";infer_eltypes=true))
# Parameter order (Total 11): λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, μE, T0_exp
df_param_sans_S = DataFrame(XLSX.readtable("../../figure_data_files/v9_files/parameters_model_1_sans_S.xlsx","Sheet1";infer_eltypes=true))

macaque_names = ["29915", "29925", "31041", "AV979",
"BA081", "BA209", "BB598", "BC094",
"BC179", "BC657", "BD536", "BD885",
"BG927", "BL669", "BO186", "BO413"]
# Progressors are placed in the top row, while controllers are placed in the remaining slots
ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]; # Order in which subplots must be placed

#=
Section 2: Define the required storage variables
=#

time_vec = 0:0.01:620; tspan = (0, 620);
N = 1e6; C = 1e6;
sol_V        = zeros(16,length(time_vec));
sol_V_sans_S = zeros(16,length(time_vec));
sol_E        = zeros(16,length(time_vec));
sol_E_sans_S = zeros(16,length(time_vec));
sol_ψ        = zeros(16,length(time_vec));
sol_ψ_sans_S = zeros(16,length(time_vec));

#=
Section 3: Solve the model for each macaque and store the values
=#

# Best-fit model

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
    sol_E[i,:] .= soln[4,:]
    sol_ψ[i,:] .= soln[4,:].*soln[5,:]
end

# Model sans S

for i = 1:1:16
    # Section 3a: Get the parameter and initial condition vectors
    curr_p    = Vector{Float64}(df_param_sans_S[i,2:12])
    init_vec  = [0.0,(10^(-2.76))/(curr_p[7]),0.0,0.0,0.0];
    curr_T0   = 10^(curr_p[11])
    curr_init = init_vec
    if i == 8 || i == 9 || i == 11 || i == 15
        curr_init[1] = curr_T0
    else
        curr_init[1] = curr_T0
        curr_init[2] = 10*curr_init[2]
    end

    # Section 3b: Solve the model
    prob = ODEProblem(model_1_sans_S_constant_k!,curr_init,tspan,curr_p)
    soln = solve(prob,saveat = 0.01,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required values and save them
    sol_V_sans_S[i,:] .= log10.(2*curr_p[7].*max.(soln[2,:], vec((1e-8)*ones(1,length(time_vec)))))
    sol_E_sans_S[i,:] .= soln[4,:]
    sol_ψ_sans_S[i,:] .= soln[4,:]
end

#=
Section 4: Plot
=#

x_lim  = [-10, 200];
y_lim  = [-2.4, 1];
x_tick = [0, 50, 100, 150, 200];
y_tick = [-2, -1, 0, 1];

fig = Figure(resolution = size_pt, fontsize=12);
order = 1; # Order counter
for i = 1:1:4
    for j = 1:1:4

        # Choosing the macaque
        idx = ordering_idx[order]; # Macaque counter

        # Setting the panel
        plot_title = macaque_names[idx];
        if plot_title == "31041" || plot_title == "AV979" || plot_title == "BB598" || plot_title == "BD885"
            ax = Axis(fig[i, j], title = plot_title, titlecolor = :red, spinewidth = 1, xticks = x_tick, yticks = y_tick,  xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
        else
            ax = Axis(fig[i, j], title = plot_title, titlecolor = :black, spinewidth = 1, xticks = x_tick, yticks = y_tick, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
        end
        hidexdecorations!(ax, ticks = false, ticklabels = false);
        hideydecorations!(ax, ticks = false, ticklabels = false);
        hidespines!(ax, :t, :r);

        # Plotting commands
        lines!(time_vec, log10.(sol_ψ[idx,:]), color = :dodgerblue1, linewidth = 2);
        lines!(time_vec, log10.(sol_ψ_sans_S[idx,:]), color = :dodgerblue1, linestyle = :dash, linewidth = 2);

        xlims!(ax, x_lim);
        ylims!(ax, y_lim);
        order += 1;
    end
end

fig
save("fig_S15a.svg", fig, px_per_unit = 6)

