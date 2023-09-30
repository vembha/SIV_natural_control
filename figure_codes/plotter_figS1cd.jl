#===========================
PLOTTER FILE FOR 
FIGURE S1C - D
===========================#

#=
Section 1: Loading the required packages and files
=#

using CairoMakie, XLSX, CSV, DataFrames
using Statistics, DifferentialEquations
CairoMakie.activate!(type="svg")
size_cm = (4.23, 4.23);
size_pt = 38 .* size_cm;
set_theme!(Theme(font = "Arial"))

#==========================
FIGURE S1C:
V VS. TIME FOR VARIOUS 
VALUES OF σ
==========================#

#=
Section 2a: Define the ODE model, the parameters, 
and solve for various values of kE
=#

function exVivo_model!(du, u, p, t)
    du[1] = -p[1]*u[1]*u[3]
    du[2] = p[2]*p[1]*u[1]*u[3] - p[3]*u[2]
    du[3] = p[4]*u[2] - p[5]*u[3]
end

β = 10^(-8); δ = 1.44; ρ = 0.36;
p = 1440; c = 0.35; f = 0.5;
ϕ = p*ρ/c; T₀ = 10^6; V₀ = 10^2.86;
kE_vec = [0, 2, 4];

u0 = [T₀; 0; V₀]; tspan = (0.0, 20.0);
p = [β, f, ρ, ϕ, δ];
prob = ODEProblem(exVivo_model!, u0, tspan, p);
sol = solve(prob);

#=
Section 3a: Define the analytical expression for V
=#

function func_V(kE, t)
    γ = kE + δ
    α = √((γ - ρ)^2 + 4*f*β*T₀*ϕ)
    V = (0.5*V₀*exp(-0.5*t*(γ + ρ))/α)*((γ - ρ + α)*(exp(-0.5*t*α)) + (ρ - γ + α)*(exp(0.5*t*α)));
    return V
end

#=
Section 4a: Solve the ODE for different values of kE
=#

V_ODE_vec = zeros(3, length(0:0.01:20));
V_fun_vec = zeros(3, length(0:0.01:20));
time_vec  = 0:0.01:20;

for i = 1:1:3
    # Obtain the function output
    current_kE = kE_vec[i]
    V_fun_vec[i, :] .= func_V.(current_kE, time_vec)

    # Solve the ODE system
    p_new = p; p_new[5] = current_kE + δ
    prob = remake(prob; p = p_new)
    sol = solve(prob, saveat = 0.01)
    V_ODE_vec[i, :] .= sol[3, :]
end

V_fun_vec = log10.(V_fun_vec);
for i = 1:1:3
    for j = 1:1:length(time_vec)
        if V_ODE_vec[i, j] < 1e-16
            V_ODE_vec[i, j] = 1e-16;
        end
    end
end
V_ODE_vec = log10.(V_ODE_vec);

#=
Section 5a: Plot
=#

fig = Figure(resolution = size_pt, fontsize=12);
ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Log viral load", xticks = [0, 10, 20]; spinewidth = 1, xtickwidth = 1, ytickwidth = 1)
hidexdecorations!(ax, ticks = false, ticklabels = false, label = true);
hideydecorations!(ax, ticks = false, ticklabels = false, label = true);
hidespines!(ax, :t, :r);

lines!(time_vec, V_ODE_vec[1,:], color = :grey, linestyle = nothing, linewidth = 2);
lines!(time_vec, V_fun_vec[1,:], color = :purple, linestyle = nothing, linewidth = 2);
lines!(time_vec, V_ODE_vec[2,:], color = :grey, linestyle = :dash, linewidth = 2);
lines!(time_vec, V_fun_vec[2,:], color = :purple, linestyle = :dash, linewidth = 2);
lines!(time_vec, V_ODE_vec[3,:], color = :grey, linestyle = :dot, linewidth = 2);
lines!(time_vec, V_fun_vec[3,:], color = :purple, linestyle = :dot, linewidth = 2);

xlims!(ax, [0, 20]); ylims!(ax, [0, 10]);
resize_to_layout!(fig);
fig
save("fig_S1c.svg", fig, px_per_unit = 6)

#==========================
FIGURE S1D:
S VS. σ FOR THE EQUATION 
AND NUMERICAL INTEGRATION 
OF THE MODEL
==========================#

#=
Section 2c: Define the parameters and analytical expression for S(t)
=#

β = 1e-8; δ = 1.44; ρ = 0.36;
p = 1440; c = 0.35; ϕ = p*ρ/c;
f = 0.5; E₀ = 10^6; T₀ = 10^6;
V₀ = 10^2.86; τₘ = 5.64;

function S(F)

    γ₄ = δ;
    γ₈ = δ + F;

    α₄ = √((γ₄ - ρ)^2 + 4*f*β*T₀*ϕ);
    α₈ = √((γ₈ - ρ)^2 + 4*f*β*T₀*ϕ);

    V₄ = (0.5*V₀*exp(-0.5*τₘ*(γ₄ + ρ))/α₄)*((γ₄ - ρ + α₄)*(exp(-0.5*τₘ*α₄)) + (ρ - γ₄ + α₄)*(exp(0.5*τₘ*α₄)));
    V₈ = (0.5*V₀*exp(-0.5*τₘ*(γ₈ + ρ))/α₈)*((γ₈ - ρ + α₈)*(exp(-0.5*τₘ*α₈)) + (ρ - γ₈ + α₈)*(exp(0.5*τₘ*α₈)));

    return log10(V₄) - log10(V₈);

end

#=
Section 3c: Define the ODE model, the parameters, 
and solve for kE = 0
=#

function exVivo_model!(du, u, p, t)
    du[1] = -p[1]*u[1]*u[3]
    du[2] = p[2]*p[1]*u[1]*u[3] - p[3]*u[2]
    du[3] = p[4]*u[2] - p[5]*u[3]
end

β = 10^(-8); δ = 1.44; ρ = 0.36;
p = 1440; c = 0.35; f = 0.5;
ϕ = p*ρ/c; T₀ = 10^6; V₀ = 10^2.86;
kE_vec = [0, 3, 6];

u0 = [T₀; 0; V₀]; tspan = (0.0, 20.0);
p = [β, f, ρ, ϕ, δ];
prob = ODEProblem(exVivo_model!, u0, tspan, p);
sol = solve(prob, saveat = 0.01); V_max_idx = findmax(sol[3, :]);
logV_withoutE = log10(V_max_idx[1]); # Maximum viral load in CD4 T cell culture
max_idx = V_max_idx[2]; # Index at which the maximum is observed

#=
Section 4c: Defining the required variables and solving the ODE system
=#

kE_vec = 0:0.1:100;
S_analytic = S.(kE_vec);
S_ODEsystem = zeros(length(kE_vec));

for i = 1:1:length(kE_vec)
    # Solving the ODE system for various kE values
    current_kE = kE_vec[i];
    p_new = p; p_new[5] = current_kE + δ
    prob = remake(prob; p = p_new);
    sol = solve(prob, saveat = 0.01);
    logV_withE = log10(sol[3, max_idx]);
    S_ODEsystem[i] = logV_withoutE - logV_withE;
end

#=
Section 5c: Plot!
=#

fig = Figure(resolution = size_pt, fontsize=12);
ax = Axis(fig[1, 1], xlabel = "F", ylabel = "S", xticks = [0, 1, 2, 3, 4, 5], yticks = [0, 1, 2, 3, 4]; spinewidth = 1, xtickwidth = 1, ytickwidth = 1)
hidexdecorations!(ax, ticks = false, ticklabels = false, label = true);
hideydecorations!(ax, ticks = false, ticklabels = false, label = true);
hidespines!(ax, :t, :r);

lines!(kE_vec, S_ODEsystem, color = :grey, linestyle = nothing, linewidth = 2);
lines!(kE_vec, S_analytic, color = :purple, linestyle = nothing, linewidth = 2);

xlims!(ax, [0, 5]); ylims!(ax, [0, 4]);
resize_to_layout!(fig);
fig
save("fig_S1d.svg", fig, px_per_unit = 6)


