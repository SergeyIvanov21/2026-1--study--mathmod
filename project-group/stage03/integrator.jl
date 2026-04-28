"""
    verlet_step!(pos, vel, forces, N, L, dt, m, ε, σ, rc)

Один шаг интегрирования методом Velocity Verlet:

  r(t+dt) = r(t) + v(t)*dt + F(t)/(2m) * dt²
  v(t+dt) = v(t) + [F(t) + F(t+dt)] / (2m) * dt

Периодические граничные условия применяются после обновления координат.
Возвращает (новые силы, потенциальная энергия на новом шаге).
"""
function verlet_step!(pos::Matrix{Float64}, vel::Matrix{Float64},
                       forces::Matrix{Float64}, N::Int, L::Float64,
                       dt::Float64, m::Float64,
                       ε::Float64, σ::Float64, rc::Float64)

    # 1. Обновление координат
    pos .+= vel .* dt .+ (forces ./ (2m)) .* dt^2

    # 2. Периодические граничные условия
    pos .= mod.(pos, L)

    # 3. Силы на новом шаге
    forces_new, U = compute_forces(pos, N, L, ε, σ, rc)

    # 4. Обновление скоростей
    vel .+= (forces .+ forces_new) ./ (2m) .* dt

    # Обновляем массив сил для следующего шага
    forces .= forces_new

    return forces_new, U
end
