#ifndef SEARCH_H
#define SEARCH_H

#include <vector>
#include <list>
#include <map>
#include <string>
#include <cstddef>

/*! \brief A trie node
 *
 *  The class implements a node of the trie structure for string matching.
 */
class Node
{
public:

    /*! \brief Node type
     *
     *  The enumeration values correspond to various pattern types.
     */
    enum Type {
        no_match = 0,   //!< a pattern corresponding to no match
        adapter = 1,    //!< an adapter fragment
        n,              //!< a gap
        polyG,          //!< a poly-G sequence
        polyC           //!< a poly-C sequence
    };

    /*! \brief The default constructor
     */
    Node() : fail(NULL), type(Type::no_match) {}
    
    /*! \brief Initialize a node with a label
     *
     *  \param[in]  label   a node label
     */
    Node(char label) :
        label(label), fail(NULL), type(Type::no_match)
    {}

    /*! \brief The node destructor
     *
     *  Deletes nodes the current node has links to.
     */
    ~Node()
    {
        for (auto it = links.begin(); it != links.end(); ++it) {
            delete *it;
        }
    }

    /*! \brief Get the next node with the specified label
     *
     *  \param[in]  c   a node label
     *  \return         the pointer to the next node with specified label
     */
    Node * next(char c)
    {
        for (size_t i = 0; i < links.size(); ++i) {
            if (links[i]->label == c)
                return links[i];
        }
        return NULL;
    }

    /*! \brief Update node statistics
     *
     *  \param[in]  t           a node type
     *  \param[in]  adapt_id    an adapter ID
     *  \param[in]  adapt_pos   an adapter position
     */
    void update_node_stats(Type t, size_t adapt_id, size_t adapt_pos)
    {
        type = t;
        adapter_id_pos.push_back(std::make_pair(adapt_id, adapt_pos));
    }

    char label;     //!< a node label
    Node * fail;    //!< a pointer to the node corresponding to matching failure
    Type type;      //!< a node type
    std::list <std::pair <size_t, size_t> > adapter_id_pos; //!< an adapter ID and position
    std::vector <Node *> links; //!< the list of links to other nodes
};

void build_trie(Node & root,
                std::vector <std::pair <std::string, Node::Type> > const & patterns,
                int errors = 0);
void add_failures(Node & root);
void go(Node * & curr, char c);
Node::Type find_match(Node * node);
Node::Type find_all_matches(Node * node, size_t pos,
                      std::map <size_t, std::vector <std::pair<size_t, size_t> > > & matches,
                      size_t errors);
Node::Type search_inexact(const std::string & text, Node * root,
                          std::vector <std::pair<std::string, Node::Type> > const & patterns, int errors);
Node::Type search_any(const std::string & text, Node * root);

#endif // SEARCH_H
