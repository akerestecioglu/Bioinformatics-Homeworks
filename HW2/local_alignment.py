

# gets the file name, mismatch, gap, match scores from the user
def get_input():
    file_name = input("Please enter the file name: ")
    mismatch = int(input("Please enter the mismatch penalty: "))
    gap = int(input("Please enter the gap penalty: "))
    match = int(input("Please enter the match score: "))
    return file_name, mismatch, gap, match


# reads the file and stores the lines in a list
def read_file(file_name):
    file = open(file_name, 'r')
    content = file.read()
    lines = content.splitlines()
    file.close()
    return lines


# creates a matrix filled with 0s
def initialize_matrix(lines):
    row_num = len(lines[0]) + 1
    col_num = len(lines[1]) + 1
    return [[0] * col_num for i in range(row_num)]


# finds the optimal local alignment by Smith-Waterman algorithm
def local_alignment(lines, row_num, col_num, gap, match, mismatch):
    mat = initialize_matrix(lines)
    trace_mat = initialize_matrix(lines)
    max_score = 0
    max_location = (0,0)

    # for every element over the list
    # note that we have one extra row and column (row 0 and col 0) filled with 0s
    for i in range(1, row_num + 1):
        for j in range(1, col_num + 1):
            options = [
                        mat[i][j - 1] + gap,
                        mat[i - 1][j] + gap,
                        mat[i - 1][j - 1] + (match if lines[0][i-1] == lines[1][j-1] else mismatch),
                        0
                      ]
            mat[i][j] = max(options)
            max_pos = options.index(max(options))
            if max_pos == 0:
                trace_mat[i][j] = (i, j-1)
            elif max_pos == 1:
                trace_mat[i][j] = (i-1, j)
            elif max_pos == 2:
                trace_mat[i][j] = (i - 1, j - 1)
            elif max_pos == 3:
                trace_mat[i][j] = (-1,-1)

            # update the maximum score if necessary
            if mat[i][j] > max_score:
                max_score = mat[i][j]
                max_location = (i, j)

    return max_score, max_location, trace_mat, mat


# traces back the score matrix and returns the local alignment
def form_alignments(max_score, max_location, trace_mat, seq1, seq2, score_mat):
    current_score = max_score
    current_loc = max_location
    align1 = ''
    align2 = ''
    while current_score > 0:
        # if we come from the diagonal cell
        if trace_mat[current_loc[0]][current_loc[1]] == (current_loc[0] - 1, current_loc[1] - 1):
            align1 = seq1[current_loc[0]-1] + align1
            align2 = seq2[current_loc[1]-1] + align2
        # if we come from the upper cell
        elif trace_mat[current_loc[0]][current_loc[1]] == (current_loc[0] - 1, current_loc[1]):
            align1 = seq1[current_loc[0]-1] + align1
            align2 = '-' + align2
        # if we come from the left cell
        elif trace_mat[current_loc[0]][current_loc[1]] == (current_loc[0], current_loc[1] - 1):
            align1 = '-' + align1
            align2 = seq2[current_loc[1]-1] + align2

        current_loc = trace_mat[current_loc[0]][current_loc[1]]
        if current_loc != (-1, -1):
            current_score = score_mat[current_loc[0]][current_loc[1]]

    return align1, align2


# writes the local alignment into the output file in the requested form
def print_alignments(align1, align2, max_score, match, mismatch, gap, lines):
    out = open('test_output5.txt', 'w')                       # output file
    mid_line = ''                                       # the middle line between sequences in the output file
    for i in range(len(align1)):
        if align1[i] == align2[i]:                      # if it is a match
            mid_line += '|'
        elif align1[i] == '-' or align2[i] == '-':      # if it is an insertion or deletion
            mid_line += ' '
        else:                                           # if it is a mismatch
            mid_line += '.'

    align1 += '\n'
    align2 += '\n'
    mid_line += '\n'
    first_line = lines[0] + '\n'
    second_line = lines[1] + '\n'
    last_line = 'Score=' + str(max_score) + ' for Match= ' + str(match) + ', mismatch = ' + str(mismatch)+ ', gap = ' + str(gap)
    out.write(first_line + second_line + '\n' + align1 + mid_line + align2 + '\n' + '\n' + last_line)
    out.close()


# main
if __name__ == '__main__':
    file_name, mismatch, gap, match = get_input()
    lines = read_file(file_name)
    max_score, max_location, trace_mat, score_mat = local_alignment(lines, len(lines[0]), len(lines[1]), gap, match, mismatch)
    align1, align2 = form_alignments(max_score, max_location, trace_mat, lines[0], lines[1], score_mat)
    print_alignments(align1, align2, max_score, match, mismatch, gap, lines)