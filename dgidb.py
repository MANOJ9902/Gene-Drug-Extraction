import requests

DGIDB_API_URL = "https://dgidb.org/api/graphql"

def query_dgidb_by_gene(query):
    """
    Query DGIdb API for drug-gene interactions based on gene name.
    """
    headers = {
        'Content-Type': 'application/json',
    }

    query_string = f'''
    {{
      genes(names: ["{query}"]) {{
        nodes {{
          interactions {{
            drug {{
              name
              conceptId
            }}
            interactionScore
            interactionTypes {{
              type
              directionality
            }}
            interactionAttributes {{
              name
              value
            }}
            publications {{
              pmid
            }}
            sources {{
              sourceDbName
            }}
          }}
        }}
      }}
    }}
    '''

    try:
        response = requests.post(DGIDB_API_URL, headers=headers, json={'query': query_string})
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while querying DGIdb: {e}")
        return {'data': {'genes': {'nodes': []}}}

def query_dgidb_by_drug(query):
    """
    Query DGIdb API for drug-gene interactions based on drug name.
    """
    headers = {
        'Content-Type': 'application/json',
    }

    query_string = f'''
    {{
      drugs(names: ["{query}"]) {{
        nodes {{
          interactions {{
            gene {{
              name
              conceptId
              longName
            }}
            interactionScore
            interactionTypes {{
              type
              directionality
            }}
            interactionAttributes {{
              name
              value
            }}
            publications {{
              pmid
            }}
            sources {{
              sourceDbName
            }}
          }}
        }}
      }}
    }}
    '''

    try:
        response = requests.post(DGIDB_API_URL, headers=headers, json={'query': query_string})
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while querying DGIdb: {e}")
        return {'data': {'drugs': {'nodes': []}}}
