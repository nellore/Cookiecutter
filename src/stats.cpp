#include "stats.h"

/*! \brief Update readc counts
 *
 *  \param[in]  type    a read type
 *  \param[in]  paired  is a read paired or not
 */
void Stats::update(ReadType type, bool paired)
{
    auto it = reads.find(type);
    if (it == reads.end()) {
        reads[type] = 1;
    } else {
        ++it->second;
    }
    ++complete;
    if (type == ReadType::ok) {
        if (paired) {
            ++pe;
        } else {
            ++se;
        }
    }
}

/*! \brief A friend function to write statistics to an output stream
 *
 *  \param[in]  out     an output stream to write read statistics to
 *  \param[in]  stats   a Stats object to write to the specified stream
 *  \return             the output stream the statistics were written to
 */
std::ostream & operator << (std::ostream & out, const Stats & stats)
{
    out << stats.filename << std::endl;
    unsigned int bad = 0;
    for (auto it = stats.reads.begin(); it != stats.reads.end(); ++it) {
        out << "\t" << get_type_name(it->first) << "\t" << it->second << std::endl;
        if (it->first != ReadType::ok) {
            bad += it->second;
        }
    }
    out << "\t" << "fraction\t" << (double)(stats.complete - bad)/stats.complete << std::endl;
    if (stats.pe) {
        out << "\t" << "se\t" << stats.se << std::endl;
        out << "\t" << "pe\t" << stats.pe << std::endl;
    }
    return out;
}
