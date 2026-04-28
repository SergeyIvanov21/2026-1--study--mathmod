using DrWatson
@quickactivate "project" 
using DifferentialEquations
using Plots

include("init.jl")
include("forces.jl")
include("integrator.jl")
include("analysis.jl")
include("visualize.jl")

# Параметры системы
const N     = 64       # Количество частиц
const L     = 10.0     # Размер ячейки моделирования
const dt    = 0.001    # Шаг по времени
const steps = 5000     # Количество шагов
const T0    = 1.0      # Начальная температура (в ед. ε/k_B)
const m     = 1.0      # Масса частицы
const rc    = 2.5      # Радиус отсечения (в ед. σ)
const ε     = 1.0      # Глубина потенциальной ямы
const σ     = 1.0      # Характерное расстояние

# Инициализация
pos   = initialize_positions(N, L)
vel   = initialize_velocities(N, T0, m)
forces, _ = compute_forces(pos, N, L, ε, σ, rc)

# Сохраняем начальные позиции для графика
pos_init = copy(pos)

# Массивы для записи результатов
E_kin = zeros(steps)
E_pot = zeros(steps)
temps = zeros(steps)

println("Запуск моделирования: N=$N частиц, $steps шагов...")

# Главный цикл
for step in 1:steps
    global forces
    forces, U = verlet_step!(pos, vel, forces, N, L, dt, m, ε, σ, rc)
    Ek          = 0.5m * sum(vel .^ 2)
    E_kin[step] = Ek
    E_pot[step] = U
    temps[step] = compute_temperature(vel, N, m)

    if step % 500 == 0
        println("Шаг $step / $steps  |  T = $(round(temps[step], digits=4))  |  E = $(round(Ek+U, digits=4))")
    end
end

println("Моделирование завершено.")

# Визуализация
plot_energy(E_kin, E_pot, steps)
plot_positions(pos_init, pos, L)
plot_velocity_distribution(vel, T0, m)
plot_temperature(temps, steps)
plot_rdf(pos, N, L, rc)

println("Графики сохранены в папку ./image/")
