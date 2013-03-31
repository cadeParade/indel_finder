from flask import Flask, render_template, request, redirect
import indel_finder_web

app = Flask(__name__)

@app.route("/")
def index():
	return render_template("index.html")

	
@app.route("/output_table", methods = ["POST"])
def make_output_page():
	b1 = request.form["b1"]
	b2 = request.form["b2"]
	gene_name = request.form["gene_name"]
	wt = request.form["wt"]
	sequences = request.form["sequences"]

	sequence_list = indel_finder_web.find_indels(b1,b2,gene_name, wt, sequences)
	table_html = indel_finder_web.outputs(sequence_list, wt)

	return render_template("table_template.html", 
						   gene_name = gene_name, data = table_html,
						   b1 = b1, b2 = b2, wt = wt, sequences = sequences)

@app.route("/csv", methods = ["GET"])
def make_csv_file():
	b1 = request.form["b1"]
	b2 = request.form["b2"]
	gene_name = request.form["gene_name"]
	wt = request.form["wt"]
	sequences = request.form["sequences"]

	csv = indel_finder_web.make_csv()

	return Response(csv, mimetype = 'text/csv')




if __name__ == "__main__":
	app.run(debug=True)
