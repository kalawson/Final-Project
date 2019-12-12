## *Background:*
PubMed is a National Library of Medicine-run database that records all publications in the greater biomedical field. PubMed assigns identities or PMIDs to all work in their database and allows users to search for keywords such as author, journal, grant number, or "disease" on their website. However, it would be useful to have an automated tracker that will mine PubMed and export the results to a useable format for later analyses.

PubMed: https://www.ncbi.nlm.nih.gov/pubmed/

## **Purpose:**
This is a tool that will automatically mine PubMed for data corresponding to a keyword of interest (i.e., grant number, author name, "cancer"). It will extract the data into a .csv file and display the information over time in bar graph format (i.e., author's publications per year or the prevelance of certain keywords over time).

The Automated PubMed Tracker tool incorporates many Python elements taught in BIOF309 including: lists, pandas data frames, loops, if/except clauses, graph generators (matplotlib), and others. We also took advantage of existing code - such as Entrez and csv converters, to produce a tool with novel PubMed search, collation, and output functions. 

This tool will be made available to other researchers interested in regularly pulling PubMed data to make bar graphs. Other potential uses include representing what journals you (or your collaborators) regularly publish in, what your publication frequency is (i.e., when are papers published), and how many last author papers you may have. Users will only need to slightly modify the keyword, area, date range in the search and tweak the matplot code to make sure your parameters are properly graphed. 

In conclusion, the Automated PubMed Tracker incorporates many elements of introductory Python coding including data visualization and manipulation. In the future, it would be interesting to parse additional paramters (such as simplifying publication date to "year") and explore alternative data visualization tools (such as Chart Studio's plotly).


## **Methods:**
Import the following modules: pandas, matplotlib.pyplot, Bio, csv, and collections. 

NLM has already developed an E-utility series called Entrez that streamlines various queries. We will be using eSearch and eSummary from the Entrez toolkit : https://dataguide.nlm.nih.gov/eutilities/utilities.html.

You will need to register with Entrez using a valid email address as part of this package.

## **Code:**
```
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez
import csv
from collections import defaultdict

Entrez.email = '------'
```

Define keyword, area, and range of years that you want to search on PubMed
```
keyword = '____'
area = '____'
start_year_range = '___'
end_year_range = '___'
```

Use eSearch to return PMIDs about your topic
```
handle = Entrez.esearch(db='pubmed', term=keyword, field = area, retmode='xml', idtype='acc', mindate=start_year_range, maxdate=end_year_range)
data = Entrez.read(handle)
```

Use eSummary to convert PMIDs to useful information
```
UID = data['IdList']
print(UID)

attributes_list = []

for ID in range(len(UID)):
    record = Entrez.esummary(db='pubmed', id=UID[ID])
    attributes_list.append(Entrez.read(record))
    
for attributes in attributes_list:
    for information in attributes:
        for key1, value1 in information.items():
            print(key1, value1)
```

Compile attributes to a .csv file
```
import csv
csv_columns = ['Item','Id','PubDate', 'EPubDate', 'Source', 'AuthorList', 'LastAuthor', 'Title', 'Volume', 'Issue', 'Pages', 'LangList', 'NlmUniqueID', 'ISSN', 'ESSN', 'PubTypeList', 'RecordStatus', 'PubStatus', 'ArticleIds', 'DOI', 'History', 'References', 'HasAbstract', 'PmcRefCount', 'FullJournalName', 'ELocationID', 'SO']
dict_data = attributes_list
csv_file = "PubMed_Output.csv"
try:
    with open(csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for info in dict_data:
            for attribute_info in info:
                writer.writerow(attribute_info)
except IOError:
    print("I/O error")
```

Format information into a citation
```
citation = defaultdict(list)

for cite in attributes_list: 
    for attribute_cite in cite:
        for key, value in attribute_cite.items():
            citation[key].append(value)

print(citation['LastAuthor'])
```

Produce bar graph from citation that counts number of references per author
```
attribute_data_plot = plt.barh(citation['LastAuthor'], citation['PmcRefCount'])
plt.ylabel = 'Number of References'
plt.xlabel = 'Author'
plt.title = 'Number of References by Author'
plt.show()
```

To see author's most cited works, search by author instead of title and use the following code:
```
author_data  = plt.barh(citation['Title'], citation['PmcRefCount'])
plt.ylabel = 'Number of References'
plt.xlabel = 'Title'
plt.title = "Most Cited Publications"
plt.show()
```

## **Authors:**
Kate Lawson (kalawson) and Taylorlyn Stephan (taylorlynstephan). Thank you to the instructors of BIOF309: John Lee, Ryan Patterson, and Sydney Hertafeld.