# particle.jl - Определение структуры частицы и операций с ней

using Random
using LinearAlgebra


# Гравитационная постоянная (реальная)
const G = 6.67430e-11

# Солнечные массы в кг
const SOLAR_MASS = 1.989e30

# Астрономическая единица в метрах
const AU = 1.496e11

# Год в секундах
const YEAR = 365.25 * 24 * 3600

# Нормированные единицы (для удобства моделирования)
# Используем: масса в массах Солнца, расстояние в АЕ, время в годах
# Тогда G_norm = G * (SOLAR_MASS * YEAR^2 / AU^3) ≈ 39.478
const G_NORM = 39.478  # В нормированных единицах (M_sun, AU, year)

# Параметры моделирования
const REPULSION_K = 1e4      # Коэффициент отталкивания
const REPULSION_A = 1.0      # Параметр отталкивания
const FRICTION_BETA = 1e3    # Коэффициент трения
const MERGE_DISTANCE = 0.05  # Расстояние для слипания

"""
    Particle
    
Структура, описывающая частицу (газопылевое облако)
"""
mutable struct Particle
    r::Vector{Float64}      # Положение [x, y, z]
    v::Vector{Float64}      # Скорость [vx, vy, vz]
    m::Float64              # Масса
    R::Float64              # Радиус
    ω::Vector{Float64}      # Угловая скорость [ωx, ωy, ωz]
    I::Float64              # Момент инерции
    active::Bool            # Активна ли частица (для слипания)
end

"""
    Particle(r, v, m, R)
    
Конструктор частицы без начального вращения
"""
function Particle(r::Vector{Float64}, v::Vector{Float64}, m::Float64, R::Float64)
    I = 0.4 * m * R^2  # I = 2/5 * m * R^2 для шара
    ω = zeros(3)
    Particle(r, v, m, R, ω, I, true)
end

"""
    Particle(r, v, m, R, ω)
    
Конструктор частицы с заданным вращением
"""
function Particle(r::Vector{Float64}, v::Vector{Float64}, m::Float64, R::Float64, ω::Vector{Float64})
    I = 0.4 * m * R^2
    Particle(r, v, m, R, ω, I, true)
end

"""
    merge_particles(p1, p2)
    
Слияние двух частиц с сохранением массы и импульса
"""
function merge_particles(p1::Particle, p2::Particle)
    @assert p1.active && p2.active "Cannot merge inactive particles"
    
    m_new = p1.m + p2.m
    R_new = cbrt(p1.R^3 + p2.R^3)
    r_new = (p1.m * p1.r + p2.m * p2.r) / m_new
    v_new = (p1.m * p1.v + p2.m * p2.v) / m_new
    
    # Сохранение момента импульса
    L1 = p1.I * p1.ω
    L2 = p2.I * p2.ω
    I_new = 0.4 * m_new * R_new^2
    ω_new = (L1 + L2) / I_new
    
    return Particle(r_new, v_new, m_new, R_new, ω_new)
end

"""
    kinetic_energy(particles)
    
Вычисление кинетической (поступательной и вращательной) энергии
"""
function kinetic_energy(particles::Vector{Particle})
    Ek_trans = 0.0
    Ek_rot = 0.0
    
    for p in particles
        if p.active
            Ek_trans += 0.5 * p.m * dot(p.v, p.v)
            Ek_rot += 0.5 * p.I * dot(p.ω, p.ω)
        end
    end
    
    return (Ek_trans, Ek_rot)
end

"""
    momentum(particles)
    
Вычисление суммарного импульса системы
"""
function momentum(particles::Vector{Particle})
    P = zeros(3)
    total_mass = 0.0
    
    for p in particles
        if p.active
            P += p.m * p.v
            total_mass += p.m
        end
    end
    
    return P, total_mass
end

"""
    zero_momentum!(particles)
    
Обнуление суммарного импульса системы
"""
function zero_momentum!(particles::Vector{Particle})
    P, total_mass = momentum(particles)
    v_cm = P / total_mass
    
    for p in particles
        if p.active
            p.v -= v_cm
        end
    end
end

"""
    count_active(particles)
    
Подсчет активных частиц
"""
function count_active(particles::Vector{Particle})
    return count(p -> p.active, particles)
end

"""
    get_active_particles(particles)
    
Получение списка активных частиц
"""
function get_active_particles(particles::Vector{Particle})
    return [p for p in particles if p.active]
end
