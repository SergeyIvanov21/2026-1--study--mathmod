using LinearAlgebra

"""
    lj_force(r_vec, ε, σ, rc)

Вычисляет силу и потенциальную энергию взаимодействия по потенциалу
Леннарда-Джонса с поправкой на сдвиг (shift) на радиусе отсечения rc.

Возвращает (F_vec, U).
"""
function lj_force(r_vec::AbstractVector, ε::Float64, σ::Float64, rc::Float64)
    r = norm(r_vec)
    if r >= rc || r < 1e-10
        return zeros(2), 0.0
    end

    sr6  = (σ / r)^6
    sr12 = sr6^2

    # Поправка на сдвиг (чтобы потенциал был непрерывен при r = rc)
    src6  = (σ / rc)^6
    src12 = src6^2
    U_shift = 4ε * (src12 - src6)

    U     = 4ε * (sr12 - sr6) - U_shift
    F_mag = 24ε / r^2 * (2 * sr12 - sr6)
    F_vec = F_mag .* r_vec

    return F_vec, U
end

"""
    compute_forces(pos, N, L, ε, σ, rc)

Вычисляет силы, действующие на все частицы, и суммарную
потенциальную энергию системы.

Использует правило минимального образа для периодических граничных условий.
"""
function compute_forces(pos::Matrix{Float64}, N::Int, L::Float64,
                         ε::Float64, σ::Float64, rc::Float64)
    forces    = zeros(N, 2)
    potential = 0.0

    for i in 1:(N - 1)
        for j in (i + 1):N
            r_vec = pos[i, :] .- pos[j, :]
            # Правило минимального образа
            r_vec .-= L .* round.(r_vec ./ L)

            F, U = lj_force(r_vec, ε, σ, rc)
            forces[i, :] .+= F
            forces[j, :] .-= F
            potential      += U
        end
    end

    return forces, potential
end
