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
    dust            //!< has low complexity according to the *DustMasker* model
};

void init_type_names(int length = 0, int polyG = 0, int dust_k = 0, int dust_cutoff = 0);

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
    Seq() : seq("") {}

    /*! \brief Read a read from a FASTQ file
     *
     *  \param[in]  fin an input stream to write the read from
     *  \return     true if a read was read, false if we reached
     *              the file end
     */
    bool read_seq(std::ifstream & fin)
    {
        std::string tmp;
        std::getline(fin, seq);
        if (!fin) {
            return false;
        }
        return true;
    }

    /*! \brief Write the read in the FASTQ format
     *
     *  \param[in]  fout    an output stream to write the read to
     */
    void write_seq(std::ofstream & fout)
    {
        fout << seq << std::endl;
    }

    /*! \brief Add the read type to its ID
     *
     *  \param[in] type the type of a read
     */
    void update_id(ReadType type)
    {
        ;
    }

    std::string seq;    //!< the read sequence
};

#endif // SEQ_H
