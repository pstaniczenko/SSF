import sys
import scipy

def compress_matrix(A):
    row_strings = []
    # create unique identifier for each row
    for i in A:
        row_strings.append(';'.join(map(str,i)))
    # find unique identifier
    unique_row_strings = list(set(row_strings))
    unique_row_strings.sort()
    # count how many of each identifier
    count = [0] * len(unique_row_strings)
    for row in row_strings:
        which_row = unique_row_strings.index(row)
        count[which_row] = count[which_row] + 1
    # retransform into scipy.array
    B = []
    for row in unique_row_strings:
        B.append(map(int, row.split(';')))
    B = scipy.array(B)
    return [B, count]

if __name__ == "__main__":
    # example:
    # python ZipData Pew2008-Q5-22429-7-Corrected.txt
    # Parameters
    adjacency_file = sys.argv[1]
    # read the data matrix
    A = scipy.genfromtxt(adjacency_file).astype(int)
    # "compress" the matrix
    B, count = compress_matrix(A)
    num_row_A = A.shape[0]
    num_col_A = A.shape[1]
    num_row_B = B.shape[0]
    # write the file
    f = open(adjacency_file + '-z-' + str(num_row_B) + '-' + str(num_col_A), 'w')
    for i in range(num_row_B):
        mystr = ' '.join(map(str, B[i])) + ' ' + str(count[i]) + '\n'
        f.write(mystr)
    f.close()
        
