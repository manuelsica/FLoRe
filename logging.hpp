// logging.hpp
// Header per il logging asincrono: definisce le variabili globali e le funzioni per gestire il log
// in un thread separato. Se viene attivata l'opzione verbose, i messaggi verranno salvati anche in un file di log.

#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <string>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>

// Utilizzando inline variables (C++17) per evitare errori di linking
inline bool logToFile = false;                   // Se true, i messaggi vengono scritti anche in un file di log
inline std::ofstream logFile;                    // File stream per il log

// Coda globale per i messaggi di log
extern std::queue<std::string> logQueue;
// Mutex per proteggere l'accesso alla coda
extern std::mutex logQueueMutex;
// Variabile condizione per notificare il thread di logging
extern std::condition_variable logQueueCV;
// Flag che indica la fine del logging
extern std::atomic<bool> loggingDone;

/*
 * async_log: Aggiunge un messaggio alla coda di log e notifica il thread dedicato.
 * Se il logging su file è abilitato, scrive anche nel file.
 */
void async_log(const std::string &msg);

/*
 * loggingThreadFunction: Funzione eseguita dal thread di logging.
 * Stampa i messaggi sulla console e, se abilitato, anche nel file di log.
 */
void loggingThreadFunction();

#endif