#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <string>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>

// Coda e variabili globali per logging asincrono
extern std::mutex logQueueMutex;
extern std::condition_variable logQueueCV;
extern std::queue<std::string> logQueue;
extern std::atomic<bool> loggingDone;

// Funzione da chiamare per accodare i messaggi di log
void async_log(const std::string &msg);

// Thread di logging principale
void loggingThreadFunction();

#endif