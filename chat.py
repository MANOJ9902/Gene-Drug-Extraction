from flask import Flask, request, render_template
from qdrant_client import QdrantClient
import openai

# Set up OpenAI API key
openai.api_key = ''
# Set up Qdrant client
client = QdrantClient(host="localhost", port=6333)

app = Flask(__name__)

def get_embedding(text):
    response = openai.Embedding.create(
        input=text,
        model="text-embedding-ada-002"
    )
    return response['data'][0]['embedding']

def get_collection_names():
    collections = client.get_collections().collections
    return [collection.name for collection in collections]

@app.route("/", methods=["GET", "POST"])
def chat():
    response = ""
    collection_names = get_collection_names()

    if request.method == "POST":
        user_input = request.form["query"]
        collection_name = request.form["collection"]

        # Convert user input to vector
        user_vector = get_embedding(user_input)

        # Retrieve relevant papers from Qdrant
        try:
            search_result = client.search(
                collection_name=collection_name,
                query_vector=user_vector,
                limit=5  # Corrected parameter
            )
        except Exception as e:
            print(f"An error occurred during the search: {e}")
            search_result = []

        # Generate response using OpenAI GPT-4
        papers_texts = " ".join([res.payload.get('abstract', '') for res in search_result])
        openai_response = openai.ChatCompletion.create(
            model="gpt-4",  # Use GPT-4 model
            messages=[
                {"role": "system", "content": "You are a helpful assistant."},
                {"role": "user", "content": f"Answer the question based on the following papers: {papers_texts}\nQuestion: {user_input}"}
            ],
            max_tokens=150,
        )
        response = openai_response.choices[0].message['content'].strip()

    return render_template("chat.html", response=response, collection_names=collection_names)

if __name__ == "__main__":
    app.run(debug=True)
