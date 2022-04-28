from scholarly import scholarly

# Schloarly Implementation (need to work on!)
def author_information(name):
    search_query = scholarly.search_author(name)
    author = scholarly.fill(next(search_query))
    [pub['bib']['title'] for pub in author['publications']]
    pub = scholarly.fill(author['publications'][0])
    return pub