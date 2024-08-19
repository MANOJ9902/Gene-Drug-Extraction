from flask import Flask, request, render_template
from biopython_query import search_protein, fetch_protein_details, fetch_papers, write_protein_info_to_file
from dgidb import query_dgidb_by_gene, query_dgidb_by_drug

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        query = request.form['query']
        query_type = request.form['query_type']

        if query_type == 'gene':
            result = query_dgidb_by_gene(query)
            pmids = []
            for interaction in result['data']['genes']['nodes']:
                for interaction_detail in interaction['interactions']:
                    for publication in interaction_detail['publications']:
                        pmids.append(str(publication['pmid']))

            papers = fetch_papers(pmids)
            return render_template('chat.html', papers=papers)
        elif query_type == 'drug':
            result = query_dgidb_by_drug(query)
            pmids = []
            for interaction in result['data']['drugs']['nodes']:
                for interaction_detail in interaction['interactions']:
                    for publication in interaction_detail['publications']:
                        pmids.append(str(publication['pmid']))

            papers = fetch_papers(pmids)
            return render_template('chat.html', papers=papers)
        elif query_type == 'protein':
            protein_ids = search_protein(query)
            protein_info = []
            for protein_id in protein_ids:
                details = fetch_protein_details(protein_id)
                if details:
                    protein_info.extend(details)
            filename = "protein_info.txt"
            write_protein_info_to_file(protein_info, filename)
            return f"Protein information written to {filename}"
        else:
            return "Invalid query type. Must be 'gene', 'drug', or 'protein'."

    return render_template('index.html')


@app.route('/chat', methods=['GET'])
def chat():
    # Example chat logic
    return render_template('chat.html')


if __name__ == "__main__":
    app.run(debug=True)
