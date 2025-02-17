import requests

def fetch_uniprot_sequence(protein_name):
    url = f'https://rest.uniprot.org/uniprotkb/search?query={protein_name}&format=fasta'
    response = requests.get(url)
    if response.status_code == 200 and response.text:
        return response.text
    return f'>No sequence found for {protein_name}\n'

def parse_proteins_from_file(input_file='genelist_unique.txt', output_file='uniprot_results.txt'):
    with open(input_file, 'r') as file:
        protein_names = [line.strip() for line in file if line.strip()]

    with open(output_file, 'w') as file:
        for protein_name in protein_names:
            sequence = fetch_uniprot_sequence(protein_name)
            file.write(sequence + '\n')
            print(sequence + '\n')
    print(f'Sequences saved to {output_file}')

# Run function with default filenames
parse_proteins_from_file()
