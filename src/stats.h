#ifndef STATS_H
#define STATS_H

#include <map>
#include <fstream>
#include <string>

#include "seq.h"

/*! \brief Statistics on processed reads
 *
 *  The class provides routines to collect and handle statistics on processed
 *  reads.
 */
class Stats
{
public:
    /*! \brief Class constructor
     *
     *  \param[in]  filename    a name of a file to write statistics to
     */
    Stats(std::string const & filename) : filename(filename), complete(0), pe(0), se(0) {}

    void update(ReadType type, bool paired = false);

    friend std::ostream & operator << (std::ostream & out, const Stats & stats);

    std::string filename;   //!< the file to write statistics to
    std::map <ReadType, unsigned int> reads;    //!< the map of read counts
    unsigned int complete;  //!< the number of processed reads
    unsigned int pe;        //!< the number of paired-end reads
    unsigned int se;        //!< the number of single-end reads
};

std::ostream & operator << (std::ostream & out, const Stats & stats);

#endif // STATS_H
