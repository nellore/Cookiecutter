#ifndef SEQ_H
#define SEQ_H

#include <string>
#include <fstream>

/*! \brief Criteria for read filtration */
enum ReadType{
    ok,             //!< passes filtration
    adapter = 1,    //!< contains an adapter sequence
    n,              //!< the number of Ns surpasses the threshold
    polyG,          //!< contains too long polyG sequence
    polyC,          //!< contains too long polyC sequence
    length,         //!< read length is too short
    dust,           //!< has low complexity according to the *DustMasker* model
    mean_quality    //!< mean quality is too small
};

void init_type_names(int length = 0, int polyG = 0, int dust_k = 0, int dust_cutoff = 0, int mean_quality = 20);

/*! \brief Get read type name from its value
 *
 *  \param[in]  type    the type of a read
 */
const std::string & get_type_name (ReadType type);

/*! \brief A sequence read with its quality information
 *
 *  This class implements routines to read reads from files in the FASTQ
 *  format.
 */
class Seq {
public:
    /*! \brief The read sequence constructor
     *
     *  It does nothing except for initializing class members
     *  (Seq::id, Seq::seq and Seq::qual) with empty strings.
     */
    Seq() : id(""), seq(""), qual("") {}

    /*! \brief Read a read from a FASTQ file
     *
     *  \param[in]  fin an input stream to write the read from
     *  \return     true if a read was read, false if we reached
     *              the file end
     */
    bool read_seq(std::ifstream & fin)
    {
        std::string tmp;
        std::getline(fin, id);
        if (!fin) {
            return false;
        }
        std::getline(fin, seq);
        std::getline(fin, tmp);
        std::getline(fin, qual);
        return true;
    }

    /*! \brief Write the read in the FASTQ format
     *
     *  \param[in]  fout    an output stream to write the read to
     */
    void write_seq(std::ofstream & fout)
    {
        fout << id << std::endl;
        fout << seq << std::endl;
        fout << '+' << std::endl;
        fout << qual << std::endl;
    }

    /*! \brief Add the read type to its ID
     *
     *  \param[in] type the type of a read
     */
    void update_id(ReadType type)
    {
        id.append(":"+get_type_name(type));
    }

    std::string id;     //!< the read ID
    std::string seq;    //!< the read sequence
    std::string qual;   //!< the read quality sequence
};

#endif // SEQ_H
