import pandas as pd

def load_pam_matrix(file_path):
    '''Load PAM scoring matrix from a file into a pandas DataFrame.'''
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f.readlines() if not line.startswith('#')]
        headers = lines[0].split()
        matrix = []
        # Skip the first element because it's the row label
        for line in lines[1:]:
            parts = line.split()
            matrix.append([int(x) for x in parts[1:]])
        pam_df = pd.DataFrame(matrix, index=headers, columns=headers)
    return pam_df

def load_fasta(file_path):
    labels, sequences = [], []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                labels.append(line.strip())
            else:
                sequences.append(line.strip())
    return labels, sequences

def init_dp(seq1, seq2, aln, gap):
    dp_scr = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    dp_dir = [[['non'] for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    if aln == 'global':
        for i in range(len(dp_scr[0])):
            dp_scr[0][i] = gap * i
            dp_dir[0][i] = ['ins']
        for i in range(len(dp_scr)):
            dp_scr[i][0] = gap * i
            dp_dir[i][0] = ['del']
        # print(*dp, sep='\n')
    return dp_scr, dp_dir

def alignment(input_path, score_path, output_path, aln, gap):
    (label1, label2), (seq1, seq2) = load_fasta(input_path)
    scores = load_pam_matrix(score_path)
    dp_scr, dp_dir = init_dp(seq1, seq2, aln, gap)

    max_score = -float('inf')
    max_poses = []
    for i in range(1, len(dp_scr)):
        for j in range(1, len(dp_scr[0])):
            c1 = seq1[j-1]
            c2 = seq2[i-1]
            score_dict = {
                'sub': dp_scr[i-1][j-1] + int(scores.loc[c2, c1]), 
                'ins': dp_scr[i][j-1] + gap,
                'del': dp_scr[i-1][j] + gap, 
                'non': 0 if aln == 'local' else -float('inf')
            }
            score = max(score_dict.values())
            dp_scr[i][j] = score
            dp_dir[i][j] = [k for k, v in score_dict.items() if v == score]
            if score > max_score:
                max_score = score
                max_poses = [(i, j)]
            elif score == max_score:
                max_poses.append((i, j))
                
            # print(dp_scr[i][j], score_dict)
            # print([k for k, v in score_dict.items() if v == dp_scr[i][j]])
    # print(*dp_scr, sep='\n')
    # print(*dp_dir, sep='\n')

    # Trace-back
    def trace_back(aln1, aln2, i, j, res):
        # print(i, j, ': ', seq1[j-1] if j > 0 else '', seq2[i-1] if i > 0 else '', dp_scr[i][j], dp_dir[i][j])
        if i <= 0 and j <= 0:
            res.append((aln1, aln2))
            return
        elif j <= 0:
            trace_back('-' + aln1, seq2[i-1] + aln2, i-1, j, res)
            return
        elif i <= 0:
            trace_back(seq1[j-1] + aln1, '-' + aln2, i, j-1, res)
            return

        if aln == 'local' and ('non' in dp_dir[i][j] or dp_scr[i][j] == 0):
            aln1 = seq1[j-1] + aln1
            aln2 = seq2[i-1] + aln2
            res.append((aln1, aln2))
            return
        if 'sub' in dp_dir[i][j]:
            trace_back(seq1[j-1] + aln1, seq2[i-1] + aln2, i-1, j-1, res)
        if 'ins' in dp_dir[i][j]:
            trace_back(seq1[j-1] + aln1, '-' + aln2, i, j-1, res)
        if 'del' in dp_dir[i][j]:
            trace_back('-' + aln1, seq2[i-1] + aln2, i-1, j, res)

    res = []
    if aln == 'global':
        i, j = len(dp_dir)-1, len(dp_dir[0])-1
        trace_back('', '', i, j, res)

    elif aln == 'local':
        for i, j in max_poses:
            trace_back('', '', i, j, res)
        max_len = max(len(a[0]) for a in res)
        res = [a for a in res if len(a[0]) == max_len]

    res.sort(key=lambda x: (x[0], x[1]))
    print(*res, sep='\n')
    with open(output_path, 'w') as f:
        for aln1, aln2 in res:
            f.write(label1 + '\n')
            f.write(aln1 + '\n')
            f.write(label2 + '\n')
            f.write(aln2 + '\n')
            if aln == 'global':
                break

if __name__ == "__main__":
    alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/my_result_global.fasta", "global", -10)
    alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/my_result_local.fasta", "local", -10)
    alignment("examples/test.fasta", "examples/pam100.txt", "examples/my_result_test.fasta", "global", -10)
    alignment("examples/input_5.fasta", "examples/pam250.txt", "examples/result5.fasta", "global", -10)
    alignment("examples/input_3.fasta", "examples/pam250.txt", "examples/result3.fasta", "local", -10)
