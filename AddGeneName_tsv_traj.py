import re

# Open the input TSV file
with open('Different_exp_patterns_across_lineages.tsv', 'r') as input_file:
    # Read the file into a list of lines
    lines = input_file.readlines()

# Open the FASTA file
with open('ICE_V1_POLISHED_MAKER2_Models.proteins.fasta', 'r') as fasta_file:
    # Read the file into a string
    fasta_str = fasta_file.read()

# Create a dictionary to store the gene names
gene_names = {}

# Use a regular expression to find the gene ID and gene name in the FASTA string
for match in re.findall(r'>Mcr_(\d+)-RA.*Name:"([^"]+)"', fasta_str):
    gene_id = 'Mcr_' + match[0]
    gene_name = match[1]
    gene_names[gene_id] = gene_name

# Add the GeneName column header to the first line of the TSV file
lines[0] = lines[0].strip() + '\tGeneName\n'

# Loop through the remaining lines of the TSV file
for line in lines[1:]:
    # Strip the line and split it into fields
    fields = line.strip().split('\t')
    # Get the gene ID from the first field
    gene_id = fields[0]
    # Get the gene name from the gene_names dictionary
    gene_name = gene_names.get(gene_id, '')
    # Add the gene name to the end of the line
    new_line = line.strip() + '\t' + gene_name + '\n'
    # Write the new line to the output file
    with open('output_file.tsv', 'a') as output_file:
        output_file.write(new_line)