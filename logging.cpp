#include "logging.hpp"
#include <iostream>

// DEFINIZIONI
bool logToFile = false;
std::ofstream logFile;

std::queue<std::string> logQueue;
std::mutex logQueueMutex;
std::condition_variable logQueueCV;
std::atomic<bool> loggingDone(false);

void async_log(const std::string &msg)
{
    {
        std::lock_guard<std::mutex> lk(logQueueMutex);
        logQueue.push(msg);
    }
    logQueueCV.notify_one();
}

void loggingThreadFunction()
{
    std::unique_lock<std::mutex> lk(logQueueMutex);
    while (true)
    {
        // aspetta nuovo messaggio o terminazione
        logQueueCV.wait(lk, []{
            return !logQueue.empty() || loggingDone.load();
        });

        // svuota la coda
        while (!logQueue.empty())
        {
            auto msg = std::move(logQueue.front());
            logQueue.pop();
            lk.unlock();

            // scrive su console
            std::cout << msg;

            // scrive su file se abilitato
            if (logToFile && logFile.is_open())
                logFile << msg;

            lk.lock();
        }

        // esci se segnato done e coda vuota
        if (loggingDone.load())
            break;
    }
}
