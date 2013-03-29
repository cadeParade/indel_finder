from flask import Flask, render_template, request, redirect
#import model
import indel_finder_web

app = Flask(__name__)

@app.route("/")
def index():
	return render_template("index.html", user_name = "GULNARA IS AWESOME", 
							another_field = "TESTING")

	
@app.route("/csv", methods = ["POST"])
def make_csv_page():
	b1 = request.form["b1"]
	b2 = request.form["b2"]
	gene_name = request.form["gene_name"]
	wt = request.form["wt"]
	sequences = request.form["sequences"]


	table_html = indel_finder_web.main(b1,b2,gene_name, wt, sequences)

	return render_template("csvtemplate.html", 
						   gene_name = gene_name, data = table_html,
						   b1 = b1, b2 = b2)



if __name__ == "__main__":
	app.run(debug=True)
