# forces.jl - Вычисление всех сил взаимодействия

using LinearAlgebra

"""
    gravity_force(p1, p2)
    
Вычисление силы гравитационного притяжения и потенциальной энергии
"""
function gravity_force(p1::Particle, p2::Particle)
    dr = p1.r - p2.r
    r = norm(dr)
    
    # Softening length для предотвращения сингулярностей
    softening = 0.5  # Увеличил для стабильности
    
    if r < 1e-10
        return zeros(3), 0.0
    end
    
    # Гравитационная сила
    F_mag = G * p1.m * p2.m / (r^2 + softening^2)
    F = F_mag * (-dr / r)
    
    # Потенциальная энергия (должна быть отрицательной!)
    potential = -G * p1.m * p2.m / sqrt(r^2 + softening^2)
    
    return (F, potential)
end

"""
    repulsion_force(p1, p2)
    
Вычисление силы отталкивания при контакте частиц
"""
function repulsion_force(p1::Particle, p2::Particle)
    dr = p1.r - p2.r
    r = norm(dr)
    R_sum = p1.R + p2.R
    
    if r >= R_sum
        return zeros(3), 0.0
    end
    
    b = r / R_sum
    F_mag = REPULSION_K * ((REPULSION_A / b)^8 - 1.0)
    F = F_mag * (dr / r)
    
    return (F, 0.0)
end

"""
    friction_force(p1, p2)
    
Вычисление силы трения и момента силы при контакте частиц
"""
function friction_force(p1::Particle, p2::Particle)
    dr = p1.r - p2.r
    r = norm(dr)
    R_sum = p1.R + p2.R
    
    if r >= R_sum
        return zeros(3), zeros(3)
    end
    
    # Относительная скорость
    W = p1.v - p2.v
    
    # Нормаль к поверхности контакта
    n = dr / r
    
    # Относительная скорость поверхностей (тангенциальная составляющая)
    # В упрощенной 2D модели для расчетов используем проекцию на плоскость
    W_perp = W - dot(W, n) * n
    
    # Для учета вращения (в 2D)
    # Вращение в плоскости XY
    ω1_effective = p1.ω[3]
    ω2_effective = p2.ω[3]
    
    # Вычисление тангенциальной скорости с учетом вращения
    # Вектор, перпендикулярный dr в плоскости XY
    n_perp = [-n[2], n[1], 0.0]
    W_perp_rot = W_perp - (ω1_effective * p1.R + ω2_effective * p2.R) * n_perp
    
    # Сила трения
    b = r / R_sum
    F_rep = REPULSION_K * ((REPULSION_A / b)^8 - 1.0)
    F_mag = FRICTION_BETA * norm(W_perp_rot) * F_rep
    
    if norm(W_perp_rot) > 0
        F_dir = -W_perp_rot / norm(W_perp_rot)
    else
        F_dir = zeros(3)
    end
    
    F = F_mag * F_dir
    
    # Момент силы
    τ = cross(dr, F)
    
    return (F, τ)
end

"""
    central_force(p, center, ω0, r0, M_center)
    
Центральная сила притяжения к неподвижной точке (Задание 1)
"""
function central_force(p::Particle, center::Vector{Float64}, ω0::Float64, r0::Float64, M_center::Float64)
    dr = center - p.r
    r = norm(dr)
    
    if r < 1e-10
        return zeros(3)
    end
    
    # Сила для поддержания круговой орбиты
    v_circ = sqrt(G * M_center / r)
    F_mag = p.m * v_circ^2 / r
    F = F_mag * (dr / r)
    
    return F
end

"""
    compute_all_forces(particles, use_gravity, use_friction)
    
Вычисление всех сил для всех частиц
"""
# forces.jl

function compute_all_forces(particles::Vector{Particle}, use_gravity::Bool, use_friction::Bool)
    n = length(particles)
    forces = [zeros(3) for _ in 1:n]
    torques = [zeros(3) for _ in 1:n]
    potential_energy = 0.0
    
    # Перебор всех пар частиц
    for i in 1:n
        if !particles[i].active
            continue
        end
        
        for j in i+1:n
            if !particles[j].active
                continue
            end
            
            # 1. Гравитационное взаимодействие
            if use_gravity
                F_grav, U_grav = gravity_force(particles[i], particles[j])
                forces[i] += F_grav
                forces[j] -= F_grav
                potential_energy += U_grav  # Потенциальная энергия пары
            end
            
            # 2. Отталкивание (всегда включено, если частицы касаются)
            F_rep, _ = repulsion_force(particles[i], particles[j])
            forces[i] += F_rep
            forces[j] -= F_rep
            
            # 3. Трение (только если включено и частицы касаются)
            if use_friction
                F_fric, τ_fric = friction_force(particles[i], particles[j])
                forces[i] += F_fric
                forces[j] -= F_fric
                torques[i] += τ_fric
                torques[j] -= τ_fric
            end
        end
    end
    
    return (forces, torques, potential_energy)
end


