// profiling.hpp
// Header per il modulo di profiling e benchmarking.
// Fornisce una classe per misurare il tempo di esecuzione e segnare eventi durante il workflow.
#ifndef PROFILING_HPP
#define PROFILING_HPP

#include <chrono>
#include <string>

// Classe Profiling per misurare il tempo e segnare eventi
class Profiling {
public:
    // Costruttore: inizializza il timer
    Profiling();

    // Avvia il profiling con una descrizione iniziale
    void start(const std::string &desc);

    // Segna un evento e stampa il tempo trascorso dall'inizio
    void mark(const std::string &label);

    // Termina il profiling e stampa il tempo totale di esecuzione
    void stop();

    // Restituisce il tempo di inizio
    std::chrono::steady_clock::time_point start_time() const;
    
private:
    std::chrono::steady_clock::time_point startTime;
};

#endif