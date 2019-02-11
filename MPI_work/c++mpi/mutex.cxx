#include <iostream>
#include <mutex>
#include <thread>
#include <string>
#include <map>

std::map<std::string, std::string> g_pages;
std::mutex g_pages_mutex;
