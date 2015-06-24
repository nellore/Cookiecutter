#include <cstddef>
#include <dirent.h>
#include <sys/stat.h>
#include "fileroutines.h"

/*! \brief Get a filename from a path
 *
 *  \param[in]  path    a file path
 *  \return             a filename from the specified path
 */
std::string basename(std::string const & path)
{
    std::string res(path);
    size_t pos = res.find_last_of('/');
    if (pos != std::string::npos) {
        res.erase(0, pos);
    }
    res = remove_extension(res);
    return res;
}

/*! \brief Remove an extension from a filename
 *
 *  \param[in]  filename    a name of a file
 *  \return                 the specified filename without its extension
 */
std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

/*! \brief Check if a directory exists, otherwise create it
 *
 *  \param[in]  dirname     a directory name
 *  \return     \p true if the specified directory exists or was created,
 *              \p false if failed to create it
 */
bool verify_directory(const std::string & dirname) {
    const char *name = dirname.c_str();
    if (opendir(name) == nullptr) {
        if (mkdir(name, 0700)) return false;
    }
    return true;
}
