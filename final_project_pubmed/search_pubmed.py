def search_pubmed(keyword, area, start_year_range, end_year_range):
    handle = Entrez.esearch(db='pubmed', term=keyword, field = area, retmode='xml', idtype='acc', mindate=start_year_range, maxdate=end_year_range)
    data = Entrez.read(handle)
    UID = data['IdList']
    attributes_list = []
    for ID in range(len(UID)):
        record = Entrez.esummary(db='pubmed', id=UID[ID])
        attributes_list.append(Entrez.read(record))
    for attributes in attributes_list:
        for information in attributes:
            for key1, value1 in information.items():
                return(key1, value1)