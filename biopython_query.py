from Bio import Entrez

# Email address required by NCBI Entrez
Entrez.email = "dolcy@cellstrat.com"


def search_protein(query):
    """
    Perform a protein search on NCBI and return the list of results.
    """
    try:
        handle = Entrez.esearch(db="protein", term=query, retmax=10, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record['IdList']
    except Exception as e:
        print(f"An error occurred while searching for proteins: {e}")
        return []


def fetch_protein_details(protein_id):
    """
    Fetch details of a protein given its ID.
    """
    try:
        handle = Entrez.efetch(db="protein", id=protein_id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record
    except Exception as e:
        print(f"An error occurred while fetching protein details: {e}")
        return None


def fetch_papers(pmids):
    """
    Fetch details of papers given a list of PMIDs.
    """
    papers = []
    try:
        handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for record in records['PubmedArticle']:
            paper_info = {}
            paper_info['pmid'] = record['MedlineCitation']['PMID']
            pmc_ids = [id_elem for id_elem in record['PubmedData']['ArticleIdList'] if
                       id_elem.attributes['IdType'] == 'pmc']
            paper_info['pmc'] = pmc_ids[0] if pmc_ids else 'N/A'
            paper_info['title'] = record['MedlineCitation']['Article']['ArticleTitle']
            paper_info['abstract'] = record['MedlineCitation']['Article']['Abstract']['AbstractText'][
                0] if 'Abstract' in record['MedlineCitation']['Article'] else 'N/A'
            paper_info['keywords'] = [kw['DescriptorName'] for kw in
                                      record['MedlineCitation'].get('MeshHeadingList', [])]
            papers.append(paper_info)
    except Exception as e:
        print(f"An error occurred while fetching papers: {e}")

    return papers


def write_protein_info_to_file(protein_info, filename):
    """
    Write protein information to a text file.
    """
    try:
        with open(filename, 'w') as f:
            for info in protein_info:
                if info:
                    f.write("Protein ID: {}\n".format(info.get('GBSeq_primary-accession', 'N/A')))
                    f.write("Description: {}\n".format(info.get('GBSeq_definition', 'N/A')))
                    f.write("Sequence Length: {}\n".format(info.get('GBSeq_length', 'N/A')))
                    f.write("-" * 50 + "\n")
        print(f"Protein information written to {filename}")
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")
