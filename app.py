#/usr/bin/env python3

from flask import Flask, render_template, request
import seqtools
from Bio.Seq import Seq
import Bio
from Bio import pairwise2
from Bio.SeqUtils import GC

#################################################################
#
# This is my fantastic Bioinformatics Web Application
#
#################################################################


# My app 
app = Flask(__name__)

@app.route('/')
def main():
    return render_template('index.html')

@app.route('/transform.html')
def transform():
    return render_template('transform.html')

@app.route('/transform_results', methods=['POST'])
def transform_results():
    input_seq = seqtools.clean(request.form['input_seq'])
    transform = request.form['transform']
    input_seq = Seq(input_seq)
    if transform == 'translate':
        output_seq = seqtools.cutseq(input_seq).translate()
    elif transform == 'transcribe':
        output_seq = input_seq.transcribe()
    output_text = seqtools.chunks(output_seq, 120)
    input_text = seqtools.chunks(input_seq, 120)
    return render_template('transform_results.html', **locals())

@app.route('/motif.html')
def motif():
    return render_template('motif.html')

@app.route('/motif_results', methods=['POST'])
def motif_results():
    input_seq = seqtools.clean(request.form['input_seq'])
    motif_seq = seqtools.clean(request.form['motif_seq'])
    input_text = seqtools.chunks(input_seq, 120)
    motif_text = seqtools.chunks(motif_seq, 120)
    if motif_seq in input_seq:
        out_text = 'The sequence DOES contain the motif'
    else:
        out_text = 'The sequence does NOT contain the motif'
    return render_template('motif_results.html', **locals())

@app.route('/complementary.html')
def complementary():
    return render_template('complementary.html')

@app.route('/complementary_results', methods=['POST'])
def complementary_results():
    input_seq = Seq(seqtools.clean(request.form['input_seq']))
    input_text=str(input_seq)
    output_text=str(input_seq.complement())

    return render_template('complementary_results.html', **locals())

@app.route('/complementary_reverse.html')
def revcomp():
    return render_template('complementary_reverse.html')

@app.route('/complementary_reverse_results', methods=['POST'])
def revcomp_results():
    input_seq = Seq(seqtools.clean(request.form['input_seq']))
    input_text= str(input_seq)
    output_text = str(input_seq.complement())
    return render_template('complementary_reverse_results.html', **locals())

@app.route('/GCcont.html')
def GCcont():
    return render_template('GCcont.html')

@app.route('/GCcont_results', methods=['POST'])
def GCcont_results():
    input_seq = Seq(seqtools.clean(request.form['input_seq']))
    input_text= str(input_seq)
    output_text = str(GC(input_seq))
    return render_template('GCcont_results.html', **locals())

@app.route('/seq_comparison.html')
def seq_comp():
    return render_template('seq_comparison.html')

@app.route('/seq_comparison_results', methods=['POST'])
def seq_comp_results():
    input_seq1 = Seq(seqtools.clean(request.form['input_seq1']))
    input_seq2 = Seq(seqtools.clean(request.form['input_seq2']))
    input_text1= str(input_seq1)
    input_text2 = str(input_seq2)
    alignments = pairwise2.align.globalxx(input_seq1, input_seq2)
    output_text=str(pairwise2.format_alignment(*alignments[0]).join(' '))
    #output_text=alignments[0]


    return render_template('seq_comparison_results.html', **locals())
# Start the app (let's keep debug=True during debugging)


app.run(debug=True)
