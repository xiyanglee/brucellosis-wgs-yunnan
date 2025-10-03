import subprocess

def run_trimal(input_file, output_file, html_out):
    subprocess.call([
        '/public/apps/trimAl-1.5.0/bin/trimal', 
        '-in', input_file, 
        '-out', output_file, 
        '-htmlout', html_out, 
        '-fasta', 
        '-gt', '0.1',  
        '-w', '1', 
        '-keepseqs', 
        '-keepheader'
    ])
