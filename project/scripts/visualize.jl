using Plots

mkpath("./image")

"""
    plot_energy(E_kin, E_pot, steps)

График кинетической, потенциальной и полной энергии.
"""
function plot_energy(E_kin::Vector{Float64}, E_pot::Vector{Float64}, steps::Int)
    t      = 1:steps
    E_tot  = E_kin .+ E_pot

    p = plot(t, E_kin, label="Кинетическая Eₖ", lw=1.5, color=:blue)
    plot!(p, t, E_pot,  label="Потенциальная U",  lw=1.5, color=:red)
    plot!(p, t, E_tot,  label="Полная E",          lw=2,   color=:green, ls=:dash)

    xlabel!(p, "Шаг")
    ylabel!(p, "Энергия (ед. ε)")
    title!(p, "Энергия системы во времени")

    savefig(p, "./image/energy.png")
    println("  → energy.png сохранён")
end

"""
    plot_positions(pos_init, pos_final, L)

Начальное и конечное расположение частиц.
"""
function plot_positions(pos_init::Matrix{Float64}, pos_final::Matrix{Float64}, L::Float64)
    p1 = scatter(pos_init[:, 1], pos_init[:, 2],
                 markersize=4, legend=false, color=:blue,
                 xlims=(0, L), ylims=(0, L),
                 title="Начальное положение", xlabel="x", ylabel="y",
                 aspect_ratio=:equal)

    p2 = scatter(pos_final[:, 1], pos_final[:, 2],
                 markersize=4, legend=false, color=:red,
                 xlims=(0, L), ylims=(0, L),
                 title="Конечное положение", xlabel="x", ylabel="y",
                 aspect_ratio=:equal)

    savefig(p1, "./image/positions_init.png")
    savefig(p2, "./image/positions_final.png")
    println("  → positions_init.png и positions_final.png сохранены")
end

"""
    plot_velocity_distribution(vel, T0, m)

Гистограмма модулей скоростей и теоретическое распределение Максвелла.
"""
function plot_velocity_distribution(vel::Matrix{Float64}, T0::Float64, m::Float64)
    speeds = sqrt.(vel[:, 1] .^ 2 .+ vel[:, 2] .^ 2)
    v_max  = maximum(speeds) * 1.1
    v_theory = range(0, v_max, length=200)
    f_theory = maxwell_boltzmann.(collect(v_theory), T0, m)

    p = histogram(speeds, normalize=:pdf, bins=20,
                  label="Моделирование", color=:skyblue, alpha=0.7)
    plot!(p, v_theory, f_theory,
          label="Максвелл (теория)", lw=2, color=:red)

    xlabel!(p, "Скорость v")
    ylabel!(p, "f(v)")
    title!(p, "Распределение скоростей частиц")

    savefig(p, "./image/velocity_dist.png")
    println("  → velocity_dist.png сохранён")
end

"""
    plot_temperature(temps, steps)

График температуры системы во времени.
"""
function plot_temperature(temps::Vector{Float64}, steps::Int)
    p = plot(1:steps, temps, lw=1.5, color=:orange, legend=false)
    hline!(p, [mean(temps)], lw=2, ls=:dash, color=:black, label="Среднее")

    xlabel!(p, "Шаг")
    ylabel!(p, "Температура T (ед. ε/k_B)")
    title!(p, "Температура системы во времени")

    savefig(p, "./image/temperature.png")
    println("  → temperature.png сохранён")
end

"""
    plot_rdf(pos, N, L, rc)

Радиальная функция распределения g(r).
"""
function plot_rdf(pos::Matrix{Float64}, N::Int, L::Float64, rc::Float64)
    r_vals, g = compute_rdf(pos, N, L; r_max=rc)

    p = plot(r_vals, g, lw=2, color=:purple, legend=false)
    hline!(p, [1.0], lw=1, ls=:dash, color=:gray)

    xlabel!(p, "r (ед. σ)")
    ylabel!(p, "g(r)")
    title!(p, "Радиальная функция распределения")

    savefig(p, "./image/rdf.png")
    println("  → rdf.png сохранён")
end
