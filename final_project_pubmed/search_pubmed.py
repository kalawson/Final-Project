from Bio import Entrez
Entrez.email = 'kalawson@vassar.edu'

def search_pubmed(keyword, area, start_year_range, end_year_range):
    """Function to search PubMed for a user defined keyword, in a user defined area and time frame.  Output is a list of up to 20 papers that fit these specifications, formated as dictionaries containing information on authors, publication date, number of citations, and other relevant attributes."""
    handle = Entrez.esearch(db='pubmed', term=keyword, field = area, retmode='xml', idtype='acc', mindate=start_year_range, maxdate=end_year_range)
    data = Entrez.read(handle)
        
    UID = data['IdList']
    
    attributes_list = []
    
    for ID in range(len(UID)):
        record = Entrez.esummary(db='pubmed', id=UID[ID])
        attributes_list.append(Entrez.read(record))

    return attributes_list
