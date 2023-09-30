#===========================
PLOTTER FILE FOR 
FIGURE S13
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

DNA_data = DataFrame(XLSX.readtable("../figure_data_files/modified_DNA_post_peak.xlsx","Sheet1";infer_eltypes=true))
DNA_data = DNA_data[!, [:animal, :time, :Y, :censoring]];
macaque_names = ["29915", "29925", "31041", "AV979",
"BA081", "BA209", "BB598", "BC094",
"BC179", "BC657", "BD536", "BD885",
"BG927", "BL669", "BO186", "BO413"]
# Progressors are placed in the top row, while controllers are placed in the remaining slots
ordering_idx = [3, 4, 7, 12, 1, 2, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]; # Order in which subplots must be placed

#=
Section 2: Estimated parameters of mono-, bi-, and tri-exponential functions
=#

#=
Section 2a: Mono-exponential parameters
D(t) = D0*exp(-r1*t)
=#

D0_mono = [3.06 3.32 4.09 4.04 3.72 3.05 3.76 3.27 3.62 3.63 3.53 3.57 3.52 3.94 3.72 3.57];
D0_mono = broadcast(^, 10, D0_mono);
r1_mono = [0.01 0.083 0.011 0.012 0.0095 0.02 0.012 0.02 0.012 0.0073 0.0081 0.011 0.0065 0.0092 0.0082 0.012];

#=
Section 2b: Bi-exponential parameters
D(t) = D0*(c1*exp(-r1*t) + (1-c1)*exp(-r2*t))
=#

D0_bi = [4 3.52 4.85 4.6 4.67 4.24 4.58 4.11 4.95 4.52 4.64 4.24 4.58 4.69 4.53 4.26];
D0_bi = broadcast(^, 10, D0_bi);
c1_bi = [0.043 0.0066 0.11 0.19 0.053 0.0067 0.089 0.019 0.014 0.09 0.049 0.14 0.073 0.12 0.1 0.11];
r1_bi = [0.0075 0.0077 0.0079 0.0082 0.0074 0.0078 0.0078 0.0078 0.0076 0.0074 0.0074 0.0079 0.0073 0.0077 0.0075 0.0078];
r2_bi = [0.1 0.11 0.1 0.1 0.1 0.11 0.1 0.1 0.096 0.11 0.1 0.1 0.1 0.1 0.1 0.11];

#=
Section 2c: Tri-exponential parameters
D(t) = D0*(c1*exp(-r1*t) + c2*exp(-r2*t) + (1-c1-c2)*exp(-r3*t))
=#

D0_tri = [3.83 3.83 4.17 4.1 4.08 3.84 3.95 3.78 4.05 3.97 3.99 3.78 4.01 4.07 3.98 3.8];
D0_tri = broadcast(^, 10, D0_tri);
c1_tri = [0.015 0.016 0.015 0.013 0.024 0.0048 0.015 0.0082 0.014 0.03 0.026 0.016 0.034 0.017 0.025 0.017];
c2_tri = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
r1_tri = [0 Inf 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
r2_tri = [0 0 Inf Inf 0 Inf 0 Inf Inf Inf Inf Inf 0 0 0 Inf];
r3_tri = [0.072 0.12 0.013 0.013 0.03 0.081 0.018 0.058 0.033 0.025 0.03 0.016 0.03 0.014 0.019 0.019];

#=
Section 3: Plot
=#

fig = Figure(resolution = size_pt, fontsize=12);
order = 1; # Order counter
for i = 1:1:4
    for j = 1:1:4

        # Choosing the macaque
        idx = ordering_idx[order]; # Macaque counter

        # Setting the panel
        plot_title = macaque_names[idx];
        if plot_title == "31041" || plot_title == "AV979" || plot_title == "BB598" || plot_title == "BD885"
            ax = Axis(fig[i, j], title = plot_title, titlecolor = :red, yticks = [0, 2.5, 5]; spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
        else
            ax = Axis(fig[i, j], title = plot_title, titlecolor = :black, yticks = [0, 2.5, 5]; spinewidth = 1, xtickwidth = 1, ytickwidth = 1, titlefont = "C:/Windows/Fonts/arialbd.ttf");
        end
        hidexdecorations!(ax, ticks = false, ticklabels = false);
        hideydecorations!(ax, ticks = false, ticklabels = false);
        hidespines!(ax, :t, :r);

        # Required dataframes
        df = filter(:animal => a -> a == macaque_names[idx], DNA_data);
        df_data = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c == 0, DNA_data);
        df_cens = filter([:animal, :censoring] => (a, c) -> a == macaque_names[idx] && c != 0, DNA_data);

        # Required functions
        func_mono(t) = log10(D0_mono[idx]*exp(-r1_mono[idx]*t));
        func_bi(t)   = log10(D0_bi[idx]*(c1_bi[idx]*exp(-r1_bi[idx]*t) + (1 - c1_bi[idx])*exp(-r2_bi[idx]*t)));
        func_tri(t)  = log10(D0_tri[idx]*(c1_tri[idx]*exp(-r1_tri[idx]*t) + c2_tri[idx]*exp(-r2_tri[idx]*t) + (1 - c1_tri[idx] - c2_tri[idx])*exp(-r3_tri[idx]*t)));

        # Plotting commands
        scatter!(df_data.time, df_data.Y, marker = :rect, strokecolor = :black, strokewidth = 1, color = :tan4, markersize = 6);
        scatter!(df_cens.time, df_cens.Y, marker = :rect, strokecolor = :black, strokewidth = 1, color = :white, markersize = 6);
        lines!(1:1:df.time[end], func_mono.(1:1:df.time[end]), color = :tan4, linestyle = :dash, linewidth = 2);
        lines!(1:1:df.time[end], func_bi.(1:1:df.time[end]), color = :tan4, linewidth = 2);
        lines!(1:1:df.time[end], func_tri.(1:1:df.time[end]), color = :tan4, linestyle = :dot, linewidth = 2);

        xlims!(ax, [0, 400]);
        ylims!(ax, [0, 5]);
        order += 1;
    end
end

fig
save("fig_S13.svg", fig, px_per_unit = 6)
