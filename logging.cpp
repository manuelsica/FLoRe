// logging.cpp
// Implementa il logging asincrono: gestisce una coda di messaggi e li scrive su console e, se abilitato, anche in un file.
#include "logging.hpp"
#include <iostream>

std::queue<std::string> logQueue;
std::mutex logQueueMutex;
std::condition_variable logQueueCV;
std::atomic<bool> loggingDone(false);

void async_log(const std::string &msg)
{
    {
        std::lock_guard<std::mutex> lock(logQueueMutex);
        logQueue.push(msg);
    }
    logQueueCV.notify_one();
    // Se il logging su file è abilitato, scrive immediatamente il messaggio anche nel file
    if (logToFile && logFile.is_open()) {
        logFile << msg << std::flush;
    }
}

void loggingThreadFunction()
{
    while (true)
    {
        std::unique_lock<std::mutex> lock(logQueueMutex);
        logQueueCV.wait(lock, [] { return !logQueue.empty() || loggingDone.load(); });
        while (!logQueue.empty())
        {
            std::string msg = logQueue.front();
            logQueue.pop();
            lock.unlock();
            std::cout << msg << std::flush;
            if (logToFile && logFile.is_open()) {
                logFile << msg << std::flush;
            }
            lock.lock();
        }
        if (loggingDone.load() && logQueue.empty())
            break;
    }
}