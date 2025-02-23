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
    # Use Levenshtein distance for better handling of insertions/deletions
    len_seq1 = len(seq1) + 1
    len_seq2 = len(seq2) + 1


    dp = [[0] * len_seq2 for _ in range(len_seq1)]

    for i in range(len_seq1):
        dp[i][0] = i
    for j in range(len_seq2):
        dp[0][j] = j

    for i in range(1, len_seq1):
        for j in range(1, len_seq2):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + cost
            )

    distance = dp[-1][-1]
    max_len = max(len(seq1), len(seq2))
    return 1 - (distance / max_len)


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



seq1 = "MLPPGSNGTAYPGQFALYQQLAQGNAVGGSAGAPPLGPSQVVTACLLTLLIIWTLLGNVLVCAAIVRSRHLRANMTNVFIVSLAVSDLFVALLVMPWKAVAEVAGYWPFGAFCDVWVAFDIMCSTASILNLCVISVDRYWAISRPFRYKRKMTQRMALVMVGLAWTLSILISFIPVQLNWHRDQAASWGGLDLPNNLANWTPWEEDFWEPDVNAENCDSSLNRTYAISSSLISFYIPVAIMIVTYTRIYRIAQVQIRRISSLERAAEHAQSCRSSAACAPDTSLRASIKKETKVLKTLSVIMGVFVCCWLPFFILNCMVPFCSGHPEGPPAGFPCVSETTFDVFVWFGWANSSLNPVIYAFNADFQKVFAQLLGCSHFCSRTPVETVNISNELISYNQDIVFHKEIAAAYIHMMPNAVTPGNREVDNDEEEGPFDRMFQIYQTSPDGDPVAESVWELDCEGEISLDKITPFTPNGFH"
file_input = "sequences.txt"
file_output = "similarity_results_P21918_DRD5_HUMAN.txt"

compare_sequences(seq1, file_input, file_output)