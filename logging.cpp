// logging.cpp
// Implementa le funzioni di logging asincrono per il tool.

#include "logging.hpp"
#include <iostream>  // per std::cout, std::flush

// Definizione delle variabili globali per il logging
std::mutex logQueueMutex;                     // Mutex per accedere alla coda
std::condition_variable logQueueCV;           // Condizione per notificare quando la coda non è vuota
std::queue<std::string> logQueue;             // Coda dei messaggi di log
std::atomic<bool> loggingDone(false);         // Flag per terminare il logging

/*
 * Funzione async_log:
 * Acquisisce il lock sulla coda, aggiunge il messaggio e notifica il thread.
 */
void async_log(const std::string &msg)
{
    {
        std::lock_guard<std::mutex> lock(logQueueMutex);
        logQueue.push(msg);  // Inserisce il messaggio nella coda
    }
    logQueueCV.notify_one(); // Notifica il thread che c'è un nuovo messaggio
}

/*
 * Funzione loggingThreadFunction:
 * Cicla finché non viene segnalato che il logging è terminato e la coda è vuota.
 */
void loggingThreadFunction()
{
    while (true)
    {
        std::unique_lock<std::mutex> lock(logQueueMutex);
        // Attende finché la coda è vuota e il flag non è impostato
        logQueueCV.wait(lock, [] { return !logQueue.empty() || loggingDone.load(); });
        // Processa tutti i messaggi presenti in coda
        while (!logQueue.empty())
        {
            std::string msg = logQueue.front();
            logQueue.pop();
            lock.unlock();  // Rilascia il lock per permettere ad altri thread di inserire messaggi
            std::cout << msg << std::flush;  // Stampa il messaggio su stdout
            lock.lock();    // Riacquisisce il lock per controllare la coda
        }
        // Se il flag di terminazione è impostato e la coda è vuota, esce dal loop
        if (loggingDone.load() && logQueue.empty())
            break;
    }
}