#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <string>
#include <queue>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <fstream>

extern bool logToFile;
extern std::ofstream logFile;

extern std::queue<std::string> logQueue;
extern std::mutex logQueueMutex;
extern std::condition_variable logQueueCV;
extern std::atomic<bool> loggingDone;

// Mette in coda un messaggio da loggare
void async_log(const std::string &msg);

// Il thread che consuma la coda e scrive su console e, se abilitato, su file
void loggingThreadFunction();

#endif // LOGGING_HPP
