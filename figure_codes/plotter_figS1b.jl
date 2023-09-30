#===========================
PLOTTER FILE FOR 
FIGURE S1B
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
CairoMakie.activate!(type="svg")
size_cm = (16.92, 9.52);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

# Read the data file
df_data = DataFrame(CSV.File("../figure_data_files/master_data_exVivo_coculture.csv"));

# Read the parameter file
# Parameter order (Total 9): β_hat_exp, δ, ρ, p, c, f, μ_exp, T0_exp, V0_exp
df_param = DataFrame(XLSX.readtable("../figure_data_files/parameters_model_exVivo_coculture.xlsx","Sheet1";infer_eltypes=true))

macaque_names = ["V401", "13523", "V20195", "13170",
"15693", "13311", "13316", "13237",
"13457", "V15729", "16834", "20973",
"20555", "20625", "20784", "20483",
"20595", "20654"];
ordering_idx = 1:1:18;

#=
Section 2: Define the model and required storage variables
=#

function exVivo_model!(du, u, p, t)
    du[1] = -(p[1]/p[7])*u[1]*u[3]
    du[2] = p[6]*(p[1]/p[7])*u[1]*u[3] - p[3]*u[2]
    du[3] = (p[3]*p[4]*p[7]/p[5])*u[2] - p[2]*u[3]
end

time_vec = 0:0.1:20; tspan = (0, 20);
sol_A = zeros(18,length(time_vec));
sol_T = zeros(18,length(time_vec));

#=
Section 3: Solve the model for each macaque and store the values
=#

for i = 1:1:18
    # Section 3a: Get the parameter and initial condition vectors
    curr_p    = Vector{Float64}(df_param[i,2:10]);
    curr_p[1] = 10^(curr_p[1]); curr_p[7] = 10^(curr_p[7]);
    init_vec  = [0.0, 0.0, 0.0];
    curr_T0   = 10^(curr_p[8]); curr_V0 = 10^(curr_p[9]);
    curr_init = init_vec; curr_init[1] = curr_T0; curr_init[3] = curr_p[7]*curr_V0;

    # Section 3b: Solve the model
    prob = ODEProblem(exVivo_model!,curr_init,tspan,curr_p[1:7])
    soln = solve(prob,saveat = 0.1,isoutofdomain=(u,p,t)->any(x->x<0,u))

    # Section 3c: Extract the required values and save them
    sol_T[i,:] .= log10.(soln[1,:])
    sol_A[i,:] .= log10.(soln[3,:])
end

#=
Section 4: Plot
=#

# Section 4a: Defining the required parameters and variables
df_plt_A = df_data;
x_lim    = [0, 15];
y_liml   = [-4, 4];
y_limr   = [0, 6.2];
x_tick   = [0, 5, 10, 15];
y_tickl  = [-4, -2, 0, 2, 4];
y_tickr  = [0, 2, 4, 6];


fig = Figure(resolution = size_pt, fontsize=9);
order = 1; # Order counter
for i = 1:1:3 # Number of rows of panels
    for j = 1:1:6 # Number of macaques per row

        # Choosing the macaque
        idx = ordering_idx[order]; # Macaque counter
        plot_title = macaque_names[idx];

        # Section 4b: Plotting for each macaque
        axl = Axis(fig[i, j], title = plot_title, titlecolor = :black, spinewidth = 1, xtickwidth = 1, ytickwidth = 1, xticks = x_tick, yticks = y_tickl, titlefont = "C:/Windows/Fonts/arialbd.ttf", yticklabelcolor = :purple);
        # axr = Axis(fig[i, j], spinewidth = 1, ytickwidth = 1, yticks = y_tickr, titlefont = "C:/Windows/Fonts/arialbd.ttf", yaxisposition = :right, yticklabelcolor = :rosybrown2);
        hidexdecorations!(axl, ticks = false, ticklabels = false); hideydecorations!(axl, ticks = false, ticklabels = false);
        # hidexdecorations!(axr); hideydecorations!(axr, ticks = false, ticklabels = false);
        hidespines!(axl, :t, :r);
        # hidespines!(axr, :t, :l);

        df_macq_A = filter(:ID => a -> a == macaque_names[idx], df_plt_A);
        sort!(df_macq_A, :time); # Sorting the dataframe in ascending order of time
        scatter!(axl, df_macq_A.time, df_macq_A.y, marker = :star4, strokecolor = :black, strokewidth = 1, color = :purple, markersize = 10);
        lines!(axl, time_vec, sol_A[idx,:], color = :purple, linewidth = 2);
        # lines!(axr, time_vec, sol_T[idx,:], color = :rosybrown2, linewidth = 2);

        xlims!(axl, x_lim);
        # xlims!(axr, x_lim);
        ylims!(axl, y_liml);
        # ylims!(axr, y_limr);
        order += 1;
        if order > 18
            break
        end
    end
end

fig
save("fig_S1b.svg", fig, px_per_unit = 6)

