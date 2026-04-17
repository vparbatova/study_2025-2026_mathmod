# quickrun.jl - Быстрый запуск всех модулей

using Plots
using Random
using LinearAlgebra

# Путь к файлам (предполагается, что все файлы в той же директории)
files = ["particle.jl", "forces.jl", "integrator.jl", "simulation.jl", 
         "visualization.jl", "tasks.jl", "main.jl"]

for file in files
    if isfile(file)
        include(file)
    else
        error("Файл $file не найден в текущей директории")
    end
end

# Запуск
main()
