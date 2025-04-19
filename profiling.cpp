// profiling.cpp
// Implementa la classe Profiling per il benchmarking del tool FLoRe utilizzando std::chrono.
#include "profiling.hpp"
#include <iostream>

Profiling::Profiling()
{
    startTime = std::chrono::steady_clock::now();
}

void Profiling::start(const std::string &desc)
{
    startTime = std::chrono::steady_clock::now();
    std::cout << "[Profiling] Inizio: " << desc << std::endl;
}

void Profiling::mark(const std::string &label)
{
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration<double>(now - startTime).count();
    std::cout << "[Profiling] " << label << ": " << elapsed << " s" << std::endl;
}

void Profiling::stop()
{
    auto now = std::chrono::steady_clock::now();
    auto total = std::chrono::duration<double>(now - startTime).count();
    std::cout << "[Profiling] Tempo totale: " << total << " s" << std::endl;
}

std::chrono::steady_clock::time_point Profiling::start_time() const
{
    return startTime;
}