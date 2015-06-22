#ifndef READROUTINES_H
#define READROUTINES_H

#include <string>
#include <vector>
#include <utility>

#include "search.h"
#include "seq.h"

double get_dust_score(std::string const & read, int k);
ReadType check_read(std::string const & read, Node * root,
                    std::vector <std::pair<std::string, Node::Type> > const & patterns,
                    unsigned int length, int dust_k, int dust_cutoff, int errors);

#endif // READROUTINES_H
