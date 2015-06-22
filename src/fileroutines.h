#ifndef FILEROUTINES_H
#define FILEROUTINES_H

#include <string>

std::string basename(std::string const & path);
std::string remove_extension(const std::string & filename);
bool verify_directory(const std::string & dirname);

#endif // FILEROUTINES_H
