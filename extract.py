import time
from Bio import Entrez
import openai
import qdrant_client
from qdrant_client.models import PointStruct, VectorParams, Distance
from xml.etree import ElementTree as ET

# NCBI Entrez email
Entrez.email = "manojreddygr62@gmail.com"

# OpenAI API key
openai.api_key = ""
# Qdrant client
client = qdrant_client.QdrantClient("http://localhost:6333")


# Function to fetch papers with retry logic
def fetch_papers(query, max_retries=3):
    try:
        search_handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        pmids = search_results.get("IdList", [])
        papers = []

        for pmid in pmids:
            retries = 0
            while retries < max_retries:
                try:
                    fetch_handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
                    fetch_results = fetch_handle.read()
                    fetch_handle.close()

                    # Parse XML data
                    root = ET.fromstring(fetch_results)

                    # Extract necessary fields
                    pmid_elem = root.find(".//PubmedArticle/MedlineCitation/PMID")
                    title_elem = root.find(".//PubmedArticle/MedlineCitation/Article/ArticleTitle")
                    abstract_elem = root.find(".//PubmedArticle/MedlineCitation/Article/Abstract/AbstractText")
                    keywords_elem = root.findall(".//PubmedArticle/MedlineCitation/KeywordList/Keyword")
                    pmc_elem = root.find(
                        ".//PubmedArticle/MedlineCitation/Article/ArticleIdList/ArticleId[@IdType='pmc']")

                    pmid_text = pmid_elem.text if pmid_elem is not None else "No PMID"
                    title_text = title_elem.text if title_elem is not None else "No Title"
                    abstract_text = abstract_elem.text if abstract_elem is not None else "No Abstract"
                    keywords_text = [kw.text for kw in keywords_elem]
                    pmc_text = pmc_elem.text if pmc_elem is not None else "No PMC ID"

                    papers.append({
                        "pmid": pmid_text,
                        "title": title_text,
                        "abstract": abstract_text,
                        "keywords": keywords_text,
                        "pmc": pmc_text
                    })
                    break  # Exit retry loop on success

                except Exception as e:
                    retries += 1
                    print(
                        f"An error occurred while fetching paper with PMID {pmid}: {e}. Retry {retries}/{max_retries}")
                    time.sleep(2)  # Wait before retrying

        return papers

    except Exception as e:
        print(f"An error occurred during the search: {e}")
        return []


# Function to create Qdrant collection and store embeddings
def create_qdrant_collection(collection_name, papers):
    try:
        # Check if collection already exists
        if collection_name in client.get_collections().collections:
            client.delete_collection(collection_name)

        # Create Qdrant collection with vector configuration
        client.create_collection(
            collection_name=collection_name,
            vectors_config={
                "size": 1536,  # Size of the embedding vector
                "distance": "Cosine"
            }
        )

        for paper in papers:
            embedding = openai.Embedding.create(
                input=paper['abstract'],
                model="text-embedding-ada-002"
            )['data'][0]['embedding']

            point = PointStruct(id=int(paper['pmid']), vector=embedding, payload=paper)
            client.upsert(collection_name=collection_name, points=[point])

        print(f"Collection {collection_name} created with {len(papers)} papers.")

    except Exception as e:
        print(f"An error occurred while creating the collection: {e}")


# Main function
def main():
    query = input("Enter gene or drug name: ")
    collection_name = query.replace(" ", "_").lower()

    papers = fetch_papers(query)

    if not papers:
        print("No papers found.")
        return

    create_qdrant_collection(collection_name, papers)


if __name__ == "__main__":
    main()
