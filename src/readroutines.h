#ifndef READROUTINES_H
#define READROUTINES_H

#include <string>
#include <vector>
#include <utility>

#include "search.h"
#include "seq.h"

double get_mean_quality(std::string const & qual);
double get_dust_score(std::string const & read, int k);
ReadType check_read(std::string const & read, std::string const & qual, Node * root,
                    std::vector <std::pair<std::string, Node::Type> > const & patterns,
                    unsigned int length, int dust_k, int dust_cutoff, int errors, int mean_quality = 0);

#endif // READROUTINES_H
