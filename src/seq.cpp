#include <map>
#include "seq.h"

/*! \brief The map object to transform read type values to their names */
std::map <ReadType, std::string> type_names;

/*! \brief Initialize the map of read types and their names
 *
 *  The goal of read type names is to indicate the reason why a read was
 *  filtered. The name contains the filter type and threshold which allows for
 *  obtaining statistics on the filtered reads. The threshold
 *  values are specified in the function arguments.
 *
 *  \param[in]  length  the minimum read length
 *  \param[in]  polyG   the minimum length of a polyG region in a read
 *  \param[in]  dust_k  this value is related to the *DustMasker* algorithm
 *  \param[in]  dust_cutoff this value is related to the *DustMasker* algorithm
 *
 *  \remark The DustMasker algorithm is described in the following paper:
 *  Morgulis, Aleksandr, E. Michael Gertz, Alejandro A. Schäffer, and Richa
 *  Agarwala. "A fast and symmetric DUST implementation to mask low-complexity
 *  DNA sequences." *Journal of Computational Biology* 13, no. 5 (2006): 1028-1040.
 */
void init_type_names(int length, int polyG, int dust_k, int dust_cutoff)
{
    type_names[ReadType::ok] = "ok";
    type_names[ReadType::adapter] = "match";
    type_names[ReadType::n] = "n";
    type_names[ReadType::polyG] = "polyG" + std::to_string(polyG);
    type_names[ReadType::polyC] = "polyC" + std::to_string(polyG);
    type_names[ReadType::length] = "length" + std::to_string(length);
    type_names[ReadType::dust] = "dust" + std::to_string(dust_k) + '_' + std::to_string(dust_cutoff);
}

/*! \brief Return a string representing a read type name
 *
 *  \param[in]  type    a read type
 *  \return the string representing the specified read type
 */
const std::string & get_type_name (ReadType type) {
    return type_names[type];
}
