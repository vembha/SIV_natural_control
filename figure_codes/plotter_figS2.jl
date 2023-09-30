#===========================
PLOTTER FILE FOR 
FIGURE S2
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
CairoMakie.activate!(type="svg")
size_cm = (18, 24);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the model file
include("model_files.jl")

# Read the data file
df_data = DataFrame(CSV.File("../figure_data_files/master_data_RNA_DNA_p27_kE0.csv"))

# Read the parameter file
# Parameter order (Total 14): λ, β_prime_exp, fD, λE⁺_exp, μI, μD, r, αE, θE, αR, θX, μE, ω_exp, T0_exp
df_param = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_2.xlsx","Sheet1";infer_eltypes=true))

macaque_names = ["29915", "29925", "31041", "AV979",
"BA081", "BA209", "BB598", "BC094",
"BC179", "BC657", "BD536", "BD885",
"BG927", "BL669", "BO186", "BO413"]
# Progressors are placed in the top row, while controllers are placed in the remaining slots
ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]; # Order in which subplots must be placed

#=
Section 2: Define the required storage variables
=#

time_vec = 0:1:620; tspan = (0, 620);
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
    curr_p    = Vector{Float64}(df_param[i,2:15])
    init_vec  = [0.0,(10^(-2.76))/(curr_p[7]),0.0,0.0,0.0];
    curr_T0   = 10^(curr_p[14])
    curr_init = init_vec
    if i == 8 || i == 9 || i == 11 || i == 15
        curr_init[1] = curr_T0
    else
        curr_init[1] = curr_T0
        curr_init[2] = 10*curr_init[2]
    end

    # Section 3b: Solve the model
    prob = ODEProblem(model_2!,curr_init,tspan,curr_p)
    soln = solve(prob,saveat = 1.0,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required values and save them
    sol_V[i,:] .= log10.(2*curr_p[7].*soln[2,:])
    sol_I[i,:] .= log10.(soln[2,:] + soln[3,:])
    sol_F[i,:] .= (soln[4,:]./N).*C.*(soln[5,:])
    sol_ψ[i,:] .= (soln[4,:]./N)

    # Section 3d: Estimate S from the predictions for F
    sol_S[i,:] .= S_func!.(sol_F[i,:])
end

#=
Section 4: Plot
=#

# Section 4a: Defining the required parameters and variables
df_plt_V = df_data[df_data.data_category .== 1, :];
df_plt_I = df_data[df_data.data_category .== 2, :];
df_plt_S = df_data[df_data.data_category .== 3, :];
df_plt_F = df_data[df_data.data_category .== 4, :];
x_lim  = [-20, 620];
y_lim  = [-0.5, 8];
x_tick = [0, 200, 400, 600];
y_tick = [-2, 0, 2, 4, 6, 8];


fig = Figure(resolution = size_pt, fontsize=12);
order = 1; # Order counter
for i = 1:1:8 # Number of rows of panels
    for j = 1:1:2 # Number of macaques per row

        # Choosing the macaque
        idx = ordering_idx[order]; # Macaque counter
        plot_title = macaque_names[idx];

        # Section 4b: Plotting for each macaque
        for k = 1:1:3 # Number of markers

            # Setting the panel
            if k == 2 # Will have the macaque's name on the top
                if plot_title == "31041" || plot_title == "AV979" || plot_title == "BB598" || plot_title == "BD885"
                    ax = Axis(fig[i, k + 3*(j - 1)], title = plot_title, titlecolor = :red, spinewidth = 1, xtickwidth = 1, ytickwidth = 1, xticks = x_tick, yticks = y_tick, titlefont = "C:/Windows/Fonts/arialbd.ttf");
                else
                    ax = Axis(fig[i, k + 3*(j - 1)], title = plot_title, titlecolor = :black, spinewidth = 1, xtickwidth = 1, ytickwidth = 1, xticks = x_tick, yticks = y_tick, titlefont = "C:/Windows/Fonts/arialbd.ttf");
                end
            else
                if plot_title == "31041" || plot_title == "AV979" || plot_title == "BB598" || plot_title == "BD885"
                    ax = Axis(fig[i, k + 3*(j - 1)], titlecolor = :red, spinewidth = 1, xtickwidth = 1, ytickwidth = 1, xticks = x_tick, yticks = y_tick, titlefont = "C:/Windows/Fonts/arialbd.ttf");
                else
                    ax = Axis(fig[i, k + 3*(j - 1)], titlecolor = :black, spinewidth = 1, xtickwidth = 1, ytickwidth = 1, xticks = x_tick, yticks = y_tick, titlefont = "C:/Windows/Fonts/arialbd.ttf");
                end
            end
            # if k == 1 && i == 8
                hidexdecorations!(ax, ticks = false, ticklabels = false);
                hideydecorations!(ax, ticks = false, ticklabels = false);
            # elseif k == 1 && i != 8
            #     hidexdecorations!(ax, ticks = false, ticklabels = true);
            #     hideydecorations!(ax, ticks = false, ticklabels = false);
            # elseif k != 1 && i == 8
            #     hidexdecorations!(ax, ticks = false, ticklabels = false);
            #     hideydecorations!(ax, ticks = false, ticklabels = true);
            # elseif k != 1 && i != 8
            #     hidexdecorations!(ax, ticks = false, ticklabels = true);
            #     hideydecorations!(ax, ticks = false, ticklabels = true);
            # end
            hidespines!(ax, :t, :r);
    
            # Modifying the dataframes and plotting
            if k == 1 # Viremia plot

                df_macq_V = filter(:animal => a -> a == macaque_names[idx], df_plt_V);
                sort!(df_macq_V, :time); # Sorting the dataframe in ascending order of time
                df_data_V = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c == 0, df_macq_V);
                df_cens_V = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c != 0, df_macq_V);

                hlines!(log10(400), color = :black, linewidth = 1, linestyle = :dash);
                scatter!(ax, df_data_V.time, df_data_V.Y, marker = :circle, strokecolor = :black, strokewidth = 1, color = :magenta, markersize = 6);
                scatter!(ax, df_cens_V.time, df_cens_V.Y, marker = :circle, strokecolor = :black, strokewidth = 1, color = :white, markersize = 6);
                lines!(ax, time_vec, sol_V[idx,:], color = :magenta, linewidth = 2);

                xlims!(ax, x_lim);
                ylims!(ax, y_lim);

            elseif k == 2 # Provirus plot

                df_macq_I = filter(:animal => a -> a == macaque_names[idx], df_plt_I);
                sort!(df_macq_I, :time); # Sorting the dataframe in ascending order of time
                df_data_I = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c == 0, df_macq_I);
                df_cens_I = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c != 0, df_macq_I);

                scatter!(ax, df_data_I.time, df_data_I.Y, marker = :rect, strokecolor = :black, strokewidth = 1, color = :tan4, markersize = 6);
                scatter!(ax, df_cens_I.time, df_cens_I.Y, marker = :rect, strokecolor = :black, strokewidth = 1, color = :white, markersize = 6);
                lines!(ax, time_vec, sol_I[idx,:], color = :tan4, linewidth = 2);

                xlims!(ax, x_lim);
                ylims!(ax, y_lim);

            elseif k == 3 # Suppressive capacity plot

                df_macq_S = filter(:animal => a -> a == macaque_names[idx], df_plt_S);
                sort!(df_macq_S, :time); # Sorting the dataframe in ascending order of time
                df_data_S = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c == 0, df_macq_S);
                df_cens_S = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c != 0, df_macq_S);

                scatter!(ax, df_data_S.time, df_data_S.Y, marker = :utriangle, strokecolor = :black, strokewidth = 1, color = :orange, markersize = 6);
                scatter!(ax, df_cens_S.time, df_cens_S.Y, marker = :utriangle, strokecolor = :black, strokewidth = 1, color = :white, markersize = 6);
                lines!(ax, time_vec, sol_S[idx,:], color = :orange, linewidth = 2);

                xlims!(ax, x_lim);
                ylims!(ax, [-0.5, 4.5]);

            end

        end
        order += 1;
    end
end

fig
save("fig_S2.svg", fig, px_per_unit = 6)

