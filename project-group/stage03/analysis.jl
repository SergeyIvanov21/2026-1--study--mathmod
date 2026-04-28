using Statistics
using LinearAlgebra

"""
    compute_temperature(vel, N, m)

Вычисляет мгновенную температуру через теорему о равнораспределении.
В 2D: <Ek> = N * k_B * T  (k_B = 1 в приведённых единицах).
"""
function compute_temperature(vel::Matrix{Float64}, N::Int, m::Float64)
    Ek = 0.5m * sum(vel .^ 2)
    return Ek / N
end

"""
    compute_rdf(pos, N, L, n_bins, r_max)

Вычисляет радиальную функцию распределения g(r).
"""
function compute_rdf(pos::Matrix{Float64}, N::Int, L::Float64;
                     n_bins::Int = 100, r_max::Float64 = L / 2)
    dr      = r_max / n_bins
    hist    = zeros(n_bins)
    density = N / L^2

    for i in 1:(N - 1)
        for j in (i + 1):N
            r_vec = pos[i, :] .- pos[j, :]
            r_vec .-= L .* round.(r_vec ./ L)
            r = norm(r_vec)
            if r < r_max
                bin = floor(Int, r / dr) + 1
                if bin <= n_bins
                    hist[bin] += 2   # учитываем оба направления (i→j и j→i)
                end
            end
        end
    end

    # Нормировка: делим на число частиц, площадь кольца и плотность
    r_vals = [(k - 0.5) * dr for k in 1:n_bins]
    for k in 1:n_bins
        area      = 2π * r_vals[k] * dr
        hist[k]  /= (N * area * density)
    end

    return r_vals, hist
end

"""
    maxwell_boltzmann(v, T, m)

Двумерное распределение Максвелла–Больцмана для модуля скорости.
"""
function maxwell_boltzmann(v::Float64, T::Float64, m::Float64)
    return (m / T) * v * exp(-m * v^2 / (2T))
end
