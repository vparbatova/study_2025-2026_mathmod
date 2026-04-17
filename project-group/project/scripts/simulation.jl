# simulation.jl - Исправленная версия

using Random

"""
    generate_keplerian_orbit(r0, M_center)
    
Генерация круговой кеплеровской орбиты
"""
# simulation.jl

const MERGE_DISTANCE = 0.8  # Было 0.5 (больше порог)

function generate_keplerian_orbit(r0::Float64, M_center::Float64)
    """
    r0 - радиус диска в АЕ
    M_center - масса центральной звезды в массах Солнца
    """
    r = r0 * sqrt(rand())
    α = 2π * rand()
    
    x = r * cos(α)
    y = r * sin(α)
    
    # Кеплеровская скорость в АЕ/год
    v_circ = sqrt(G_NORM * M_center / r)
    vx = -v_circ * sin(α)
    vy = v_circ * cos(α)
    
    println("    r = $(round(r, digits=2)) AU, v = $(round(v_circ, digits=2)) AU/год")
    
    return ([x, y, 0.0], [vx, vy, 0.0], r)
end

"""
    generate_disk_particles_2d(N, r0, M_center)
    
Генерация частиц в 2D диске с правильными скоростями
"""
function generate_disk_particles_2d(N::Int, r0::Float64, M_center::Float64)
    particles = Particle[]
    
    # Суммарная масса частиц
    total_mass = 1.0
    
    println("  Генерация частиц с кеплеровскими скоростями...")
    
    for i in 1:N
        r, v, radius = generate_keplerian_orbit(r0, M_center)
        
        # Масса частицы
        m = total_mass / N
        R = cbrt(m) * 0.5  # Увеличил радиус для лучшей видимости
        
        println("    Частица $i: r=$(round(radius, digits=2)), v=$(round(norm(v), digits=4))")
        push!(particles, Particle(r, v, m, R))
    end
    
    # Вывод суммарного импульса ДО обнуления
    P_before, _ = momentum(particles)
    println("  Суммарный импульс ДО обнуления: $(round(norm(P_before), digits=6))")
    
    # Обнуление суммарного импульса
    zero_momentum!(particles)
    
    P_after, _ = momentum(particles)
    println("  Суммарный импульс ПОСЛЕ обнуления: $(round(norm(P_after), digits=6))")
    
    return particles
end

"""
    generate_disk_particles_3d(N, r0, M_center, disk_thickness)
    
Генерация частиц в 3D диске
"""
function generate_disk_particles_3d(N::Int, r0::Float64, M_center::Float64, disk_thickness::Float64=0.1)
    particles = Particle[]
    total_mass = 1.0
    
    println("  Генерация 3D частиц...")
    
    for i in 1:N
        # Радиус в плоскости диска
        r_xy = r0 * sqrt(rand())
        α = 2π * rand()
        x = r_xy * cos(α)
        y = r_xy * sin(α)
        
        # Вертикальная координата
        z = (rand() - 0.5) * disk_thickness * r0
        
        # Расстояние от центра
        r = sqrt(x^2 + y^2 + z^2)
        
        # Кеплеровская скорость
        v_circ = sqrt(G_NORM * total_mass / r)
        vx = -v_circ * sin(α)
        vy = v_circ * cos(α)
        vz = 0.0
        
        m = total_mass / N
        R = cbrt(m) * 0.5
        
        push!(particles, Particle([x, y, z], [vx, vy, vz], m, R))
    end
    
    zero_momentum!(particles)
    return particles
end

"""
    generate_two_mass_particles(N, r0, light_fraction, mass_ratio)
    
Генерация частиц двух сортов
"""
function generate_two_mass_particles(N::Int, r0::Float64, light_fraction::Float64, mass_ratio::Float64)
    particles = Particle[]
    
    n_light = Int(round(N * light_fraction))
    n_heavy = N - n_light
    
    # Распределение масс (тяжелые частицы массивнее)
    total_mass_heavy = 0.7 * N / 30  # Масштабирование
    total_mass_light = 0.3 * N / 30
    
    m_heavy = total_mass_heavy / max(n_heavy, 1)
    m_light = total_mass_light / max(n_light, 1)
    
    println("  Тяжелые частицы: $n_heavy (m=$(round(m_heavy, digits=4)))")
    println("  Легкие частицы: $n_light (m=$(round(m_light, digits=4)))")
    
    # Тяжелые частицы
    for i in 1:n_heavy
        r, v, _ = generate_keplerian_orbit(r0, 1.0)
        R = cbrt(m_heavy) * 0.6
        push!(particles, Particle(r, v, m_heavy, R))
    end
    
    # Легкие частицы
    for i in 1:n_light
        r, v, _ = generate_keplerian_orbit(r0, 1.0)
        R = cbrt(m_light) * 0.4
        push!(particles, Particle(r, v, m_light, R))
    end
    
    # Перемешивание
    particles = shuffle(particles)
    zero_momentum!(particles)
    
    return particles
end
