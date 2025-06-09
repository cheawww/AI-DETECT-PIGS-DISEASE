from flask import Flask, render_template, request

app = Flask(__name__)

def generate_grna(sequence, pam="NGG"):
    """
    Function to generate gRNAs from the given sequence.
    Looks for PAM sites (default: NGG for SpCas9).
    """
    sequence = sequence.upper()
    grna_list = []
    
    for i in range(len(sequence) - 23):
        candidate = sequence[i : i + 23]
        if candidate[-2:] == "GG":  # PAM check for SpCas9
            grna_list.append({"grna": candidate, "start": i, "end": i + 23})
    
    return grna_list

@app.route("/", methods=["GET", "POST"])
def hub():
    return render_template("hub.html")
    


@app.route("/analyze", methods=["POST"])
def analyze():
    # Extract form inputs
    genome = request.form.get("genome")
    region = request.form.get("region")
    gene = request.form.get("gene")
    chrom = request.form.get("chrom")
    sequence = request.form.get("sequence")
    email = request.form.get("email")
    jobname = request.form.get("jobname")

    # Generate gRNAs from the sequence
    grna_list = generate_grna(sequence)

    # Render hub2.html with gRNA data
    return render_template(
        "hub2.html",
        genome=genome,
        region=region,
        gene=gene,
        chrom=chrom,
        sequence=sequence,
        grna_list=grna_list,
        enumerate=enumerate,  # Add enumerate to the context
    )



if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=5500)
