# %%
from Bio import Entrez
Entrez.email = 'learnbiopython@gmail.com'

# %%
term1 = ['Boss saurus', 'Antelope cervicapra', 'Gazella bbennettii', 
         'Boselaphus tragocamelus', 'Canis lupus', 'Panthera leo', 
         'Elephas laximus', 'Equus africanus', 'Panthera pardus', 
         'Cervus canadensis', 'Pavo cristatus', 'Grus leucogeranus', 
         'Vulpes vulpes', 'Rhinoceros unicornis', 'Panthera Tigris', 
         'Crocodylus palustris', 'Gavialis gangeticus', 'Equus caballus', 
         'Equus quagga', 'Babalus bubalis', 'Sus scrofa', 'Camelus dromedaries', 
         'Giraffa camelopardalis ', 'Hemidactylus flaviviridis', 'Hippopotamus amphibius', 
         'Macaca mulatta', 'Canis lupus', 'Felis domesticus', 'Acinonyx jubatus', 
         'Rattus rattus', 'Mus musculus', 'Oryctolagus cuniculus', 'Bubo virginianus', 
         'Passer domesticus', 'Corvus splendens', 'Acridotheres tristis', 'Psittacula eupatria', 
         'Molpastes caferr', 'Eudynamis scolopaccus', 'Columba livia', 'Naja naja', 'Ophiophagus hannah', 
         'Hydrophiinae ', 'Python molurus', 'Ptyas mucosa']

term2 = [
    'Bos gaurus', 'Gazella bennettii', 'Boselaphus tragocamelus',
    'Canis lupus', 'Elephas maximus',  
    'Cervus canadensis', 'Pavo cristatus', 'Grus leucogeranus', 'Vulpes vulpes',
    'Crocodylus palustris', 'Gavialis gangeticus',
    'Equus caballus', 'Equus quagga', 'Sus scrofa', 'Camelus dromedaries',
    'Canis lupus', 'Felis domesticus','Rattus rattus',
    'Mus musculus', 'Oryctolagus cuniculus'
]
term3 = ['Tursiops truncatus']


# %%

def checkNames(term):
    handle = Entrez.espell(db='pmc', term=term)
    record = Entrez.read(handle)
    handle.close()
    correctTerm = record.get('CorrectedQuery', term)
    words = correctTerm.split()
    correctWord = [word.capitalize() if i == 0 else word.lower() for i, word in enumerate(words)]
    return " ".join(correctWord)

print("_" * 82)
print(f" | {'Given Names':^30} | {'Binomial Nomenclature':^30} | {'Status':^10} |")
print("_" * 82)

for name in term1:
    correctTerm = checkNames(name)
    status = "True" if name == correctTerm else "False"
    correctDisplay = correctTerm if status == "False" else name
    print(f" | {name:^30} | {correctDisplay:^30} | {status:^10} |")
print("_" * 82)


# %%

def commonNames(terms):

    handle = Entrez.esearch(db='taxonomy', term=terms)
    record = Entrez.read(handle)
    handle.close()

    if record['IdList']:
        tax_id = record['IdList'][0]
        handle = Entrez.efetch(db='taxonomy', id=tax_id)
        tax_record = Entrez.read(handle)
        handle.close()

        if 'CommonName' in tax_record[0]['OtherNames']:
            common_names = tax_record[0]['OtherNames']['CommonName']
            if isinstance(common_names, list):
                return ", ".join(common_names)
            else:
                return common_names
    else:
        return "N/A"

print("_" * 110)
print(f" | {'Nama Ilmiah':^31} | {'Nama Umum':^71} |")
print("_" * 110)

for name in term2:
    common_name = commonNames(name)
    print(f" | {name:^31} | {common_name:^71} |")


# %%
print('\n')
print("Common Names : ",commonNames(term3))
def fetch_taxonomy_info(terms):
    try:
        handle = Entrez.esearch(db='taxonomy', term=terms)
        record = Entrez.read(handle)
        handle.close()
        if record['IdList']:
            tax_id = record['IdList'][0]
            handle = Entrez.efetch(db='taxonomy', id=tax_id)
            tax_record = Entrez.read(handle)
            handle.close()
            return tax_record[0]
    except Exception as e:
        print(f"Error fetching taxonomy information for {terms}: {e}")
        return None

def fetch_genomic_publications(terms):
    try:
        record = Entrez.read(Entrez.esearch(db='pubmed', term=f"{terms} genome", retmax=5))
        if record['IdList']:
            publications = []
            for pub_id in record['IdList']:
                pub_summary = Entrez.read(Entrez.esummary(db='pubmed', id=pub_id))
                publications.append({
                    'Title': pub_summary[0]['Title'],
                    'FullJournalName': pub_summary[0]['FullJournalName'],
                    'DOI': pub_summary[0].get('DOI', 'N/A')
                })
            return publications
    except Exception as e:
        print(f"Error fetching genomic publications for {terms}: {e}")
        return []

def fetch_genome_database_size(terms):
    try:
        record = Entrez.read(Entrez.esearch(db='nucleotide', term=f"{terms} genome", retmax=1))
        return record['Count']
    except Exception as e:
        print(f"Error fetching genome database size for {terms}: {e}")
        return '0'

def fetch_genbank_info(terms):
    try:
        record = Entrez.read(Entrez.esearch(db='nucleotide', term=f'{terms}[Organism] AND refseq[filter]', retmax=1, idtype='acc'))
        if record['IdList']:
            ID = record['IdList'][0]
            if 'NM_' in ID:
                fetch = Entrez.efetch(db='nucleotide', id=ID, rettype='fasta', retmode='text')
                read_fetch = fetch.read()
                return read_fetch
    except Exception as e:
        print(f"Error fetching GenBank information for {terms}: {e}")
        return None

def process_termss(termss):
    for name in termss:
        tax_info = fetch_taxonomy_info(name)
        if tax_info:
            print(f"Scientific Name: {tax_info['ScientificName']}")
            print(f"Taxonomy Rank: {tax_info['Rank']}")

            lineage_list = tax_info['Lineage'].split("; ")
            print("\nLineage:")
            print("_" * 37)
            print(f" | {'Level':^5} | {'Taxon':^24} |")
            print("_" * 37)
            for level, taxon in enumerate(lineage_list, 1):
                print(f" | {level:^5} | {taxon:^25}|")
            print("_" * 37)

            print("\nRelated Genomic Publications:")
            publications = fetch_genomic_publications(name)
            for idx, pub in enumerate(publications, 1):
                print(f"{idx}. Title: {pub['Title']}")
                print(f"   Full Journal Name: {pub['FullJournalName']}")
                print(f"   DOI: {pub['DOI']}\n")

            genome_size = fetch_genome_database_size(name)
            print(f"Genome Database Size: {genome_size} entries")

            genbank_info = fetch_genbank_info(name)
            if genbank_info:
                print("\nGenBank Information:")
                print(genbank_info)
        else:
            print(f"No taxonomy information found for {name}.")

# Example list of scientific names to process
process_termss(term3)

# Fetch records from PMC
try:
    record = Entrez.read(Entrez.esearch(db='pmc', term='Felis domesticus', retmax=10))
    print(type(record))
    print(record.keys())
    for key in record.keys():
        print(key, ':', record[key])
except Exception as e:
    print(f"Error fetching PMC records: {e}")



