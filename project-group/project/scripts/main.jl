# main.jl - Главный файл программы

using Plots
using Random
using LinearAlgebra

# Подключение модулей
include("particle.jl")
include("forces.jl")
include("integrator.jl")
include("simulation.jl")
include("visualization.jl")

function main()
    println("="^60)
    println("Моделирование образования планетной системы")
    println("="^60)
    println("\nВыберите задание:")
    println("  1. Движение N точек с притяжением к центру")
    println("  2. Гравитационное взаимодействие + отталкивание")
    println("  3. Собственное вращение частиц")
    println("  4. 3D модель без трения")
    println("  5. 3D модель со слипанием")
    println("  6. Полная модель с трением")
    println("  7. Частицы двух сортов (разные массы)")
    println("  0. Все задания последовательно")
    print("\nВаш выбор: ")
    
    choice = parse(Int, readline())
    
    if choice == 0
        run_all_tasks()
    elseif choice == 1
        run_task1()
    elseif choice == 2
        run_task2()
    elseif choice == 3
        run_task3()
    elseif choice == 4
        run_task4()
    elseif choice == 5
        run_task5()
    elseif choice == 6
        run_task6()
    elseif choice == 7
        run_task7()
    else
        println("Неверный выбор!")
    end
end

function run_all_tasks()
    println("\n>>> Запуск всех заданий по порядку <<<\n")
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 1")
    println("-"^40)
    run_task1()
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 2")
    println("-"^40)
    run_task2()
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 3")
    println("-"^40)
    run_task3()
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 4")
    println("-"^40)
    run_task4()
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 5")
    println("-"^40)
    run_task5()
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 6")
    println("-"^40)
    run_task6()
    
    println("\n" * "-"^40)
    println("ЗАДАНИЕ 7")
    println("-"^40)
    run_task7()
    
    println("\n" * "="^60)
    println("Все задания выполнены!")
    println("="^60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
