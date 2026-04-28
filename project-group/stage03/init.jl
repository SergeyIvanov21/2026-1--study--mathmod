using Statistics

"""
    initialize_positions(N, L)

Размещает N частиц в узлах квадратной решётки в ячейке размером L×L.
"""
function initialize_positions(N::Int, L::Float64)
    n_side   = ceil(Int, sqrt(N))
    spacing  = L / n_side
    positions = zeros(N, 2)
    for i in 1:N
        ix = (i - 1) % n_side
        iy = (i - 1) ÷ n_side
        positions[i, 1] = (ix + 0.5) * spacing
        positions[i, 2] = (iy + 0.5) * spacing
    end
    return positions
end

"""
    initialize_velocities(N, T0, m)

Генерирует случайные скорости из распределения Максвелла.
Суммарный импульс системы обнуляется.
"""
function initialize_velocities(N::Int, T0::Float64, m::Float64)
    vel = randn(N, 2)
    # Обнуление суммарного импульса
    vel .-= mean(vel, dims=1)
    # Масштабирование на заданную температуру:
    # <m*v^2/2> = k_B*T  (2D: два степени свободы на частицу)
    Ek_current = 0.5m * sum(vel .^ 2)
    T_current  = Ek_current / N
    vel      .*= sqrt(T0 / T_current)
    return vel
end
