import csv

# Open the CSV file

file='DEG_salt_vs_control.csv'

with open(file, 'r') as csv_file:
    reader = csv.reader(csv_file)
    rows = list(reader)

    # Modify the header row to add the new column
    header = rows[0]
    header.append('Description')

    # Iterate over the rows (excluding the header)
    for row in rows[1:]:
        gene_id = row[0]

        # Open the FASTA file
        with open('ICE_V1_POLISHED_MAKER2_Models.proteins.fasta', 'r') as fasta_file:
            fasta_lines = fasta_file.readlines()

            # Search for the corresponding gene_id in the FASTA file
            for i in range(len(fasta_lines)):
                if fasta_lines[i].startswith('>' + gene_id):
                    # Extract the description text
                    description = fasta_lines[i].split('Similar to', 1)[-1].strip()

                    # Exclude text after and including "AED" in the description
                    description = description.split('AED', 1)[0].strip()

                    # Add the description to the current row
                    row.append(description)
                    break

# Write the modified CSV file
with open(file, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(rows)
