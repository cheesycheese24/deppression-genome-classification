import re


def read_fasta(file_input):
    sequences = []
    with open(file_input, 'r') as f:
        header = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>sp'):
                if header and sequence:
                    sequences.append((header, ''.join(sequence)))
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header and sequence:
            sequences.append((header, ''.join(sequence)))
    return sequences


def calculate_similarity_score(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / max(len(seq1), len(seq2))


def compare_sequences(seq1, file_input, file_output):
    sequences = read_fasta(file_input)
    results = []

    for header, seq2 in sequences:
        score = calculate_similarity_score(seq1, seq2)
        results.append((score, header))

    results.sort(reverse=True, key=lambda x: x[0])

    with open(file_output, 'w') as f:
        for score, header in results:
            f.write(f"{header}\nScore: {score:.4f}\n\n")


# Example usage
seq1 = "MAHRGPSRASKGPGPTARAPSPGAPPPPRSPRSRPLLLLLLLLGACGAAGRSPEPGRLGPHAQLTRVPRSPPAGRAEPGGGEDRQARGTEPGAPGPSPGPAPGPGEDGAPAAGYRRWERAAPLAGVASRAQVSLISTSFVLKGDATHNQAMVHWTGENSSVILILTKYYHADMGKVLESSLWRSSDFGTSYTKLTLQPGVTTVIDNFYICPTNKRKVILVSSSLSDRDQSLFLSADEGATFQKQPIPFFVETLIFHPKEEDKVLAYTKESKLYVSSDLGKKWTLLQERVTKDHVFWSVSGVDADPDLVHVEAQDLGGDFRYVTCAIHNCSEKMLTAPFAGPIDHGSLTVQDDYIFFKATSANQTKYYVSYRRNEFVLMKLPKYALPKDLQIISTDESQVFVAVQEWYQMDTYNLYQSDPRGVRYALVLQDVRSSRQAEESVLIDILEVRGVKGVFLANQKIDGKVMTLITYNKGRDWDYLRPPSMDMNGKPTNCKPPDCHLHLHLRWADNPYVSGTVHTKDTAPGLIMGAGNLGSQLVEYKEEMYITSDCGHTWRQVFEEEHHILYLDHGGVIVAIKDTSIPLKILKFSVDEGLTWSTHNFTSTSVFVDGLLSEPGDETLVMTVFGHISFRSDWELVKVDFRPSFSRQCGEEDYSSWELSNLQGDRCIMGQQRSFRKRKSTSWCIKGRSFTSALTSRVCECRDSDFLCDYGFERSSSSESSTNKCSANFWFNPLSPPDDCALGQTYTSSLGYRKVVSNVCEGGVDMQQSQVQLQCPLTPPRGLQVSIQGEAVAVRPGEDVLFVVRQEQGDVLTTKYQVDLGDGFKAMYVNLTLTGEPIRHRYESPGIYRVSVRAENTAGHDEAVLFVQVNSPLQALYLEVVPVIGLNQEVNLTAVLLPLNPNLTVFYWWIGHSLQPLLSLDNSVTTRFSDTGDVRVTVQAACGNSVLQDSRVLRVLDQFQVMPLQFSKELDAYNPNTPEWREDVGLVVTRLLSKETSVPQELLVTVVKPGLPTLADLYVLLPPPRPTRKRSLSSDKRLAAIQQVLNAQKISFLLRGGVRVLVALRDTGTGAEQLGGGGGYWAVVVLFVIGLFAAGAFILYKFKRKRPGRTVYAQMHNEKEQEMTSPVSHSEDVQGAVQGNHSGVVLSINSREMHSYLVS"
file_input = "sequences.txt"
file_output = "similarity_results.txt"

compare_sequences(seq1, file_input, file_output)
