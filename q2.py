import numpy as np

STOP, LEFT_TRACEBACK, UP_TRACEBACK, DIAG_TRACEBACK = 0, 1, 2, 3

def create_score_matrix(seq1, seq2, match = 2, mismatch = -2, gap = -3):
    """
    Cria as matrizes para representar os alinhamentos locais 
    de duas sequências usando o algoritmo de Smith-Waterman.

    Args:
        seq1 (str): A primeira sequência.
        seq2 (str): A segunda sequência.
        match_score (int): Pontuação para um match.
        mismatch_score (int): Pontuação para um mismatch.
        gap_penalty (int): Penalidade para um gap.

    Returns:
        tuple: Matriz de pontuação, matriz de direções,
               melhor valor local e posição do melhor valor
    """
    # Dimensões M e N das matrizes
    m = len(seq1) + 1
    n = len(seq2) + 1

    # Cria as matrizes para pontuação e direções
    score_matrix = np.zeros((m, n), dtype=np.int32)
    traceback_matrix = np.zeros((m, n), dtype=np.int32)

    # Inicia para encontrar os valores do melhor alinhamento depois
    max_value = 0
    max_pos = None

    for i in range(1, m):
        for j in range(1, n):
            # Verifica se teve match ou mismatch das sequências
            if seq1[i - 1] == seq2[j - 1]:
                match_value = match
            else:
                match_value = mismatch

            # Valores possiveis
            diag = score_matrix[i-1][j-1] + match_value
            up   = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap

            # Seleciona o melhor valor e marca de onde veio na matriz de direções
            if diag >= max(up, left, 0):
                score_matrix[i][j] = diag
                traceback_matrix[i][j] = DIAG_TRACEBACK
            elif up >= max(diag, left, 0):
                score_matrix[i][j] = up
                traceback_matrix[i][j] = UP_TRACEBACK
            elif left >= max(diag, up, 0):
                score_matrix[i][j] = left
                traceback_matrix[i][j] = LEFT_TRACEBACK
            else:
                score_matrix[i][j] = 0
                traceback_matrix[i][j] = STOP

            # Encontra o melhor valor e sua posição
            if score_matrix[i][j] > max_value:
                max_value = score_matrix[i][j]
                max_pos = (i, j)

    return score_matrix, traceback_matrix, max_pos, max_value


def traceback(seq1, seq2, score_matrix, traceback_matrix, start_pos):
    """
    Retorna os trechos com melhor alinhamento local.

    Args:
        seq1 (str): A primeira sequência.
        seq2 (str): A segunda sequência.
        score_matrix (matriz de ints): Matriz das pontuações
        traceback_matrix (matriz de ints): Matriz das setas para o traceback
        start_pos (int): Posição onde inicia o traceback

    Returns:
        tuple: As sequências que tiveram o melhor alinhamento local
    """
    # Inicia as variáveis para salver as sequências
    aligned1 = []
    aligned2 = []
    i, j = start_pos

    # Vai lendo e salvando a partir da posição dada até o inicio da sequência e salvando as strings
    while traceback_matrix[i][j] != STOP and score_matrix[i][j] > 0:
        if traceback_matrix[i][j] == DIAG_TRACEBACK:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == UP_TRACEBACK:
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            i -= 1
        elif traceback_matrix[i][j] == LEFT_TRACEBACK:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2))


def smith_waterman(seq1, seq2,match = 2, mismatch = -2, gap = -3):
    """
    Realiza os alinhamentos locais de duas sequências usando o algoritmo de Smith-Waterman.

    Args:
        seq1 (str): A primeira sequência.
        seq2 (str): A segunda sequência.
        match_score (int): Pontuação para um match.
        mismatch_score (int): Pontuação para um mismatch.
        gap_penalty (int): Penalidade para um gap.

    Returns:
        tuple: Matriz de pontuação, matriz de direções,
               melhor valor local e posição do melhor valor
    """
    score_matrix, traceback_matrix, start_pos, max_score = create_score_matrix(seq1, seq2, match, mismatch, gap)
    aligned1, aligned2 = traceback(seq1, seq2, score_matrix, traceback_matrix, start_pos)
    return int(max_score), aligned1, aligned2


def main():
    # Dados para as sequências de entrada
    seq_hs = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKY"
    seq_chrysocyon = "VLSPADKTNIKSTWDKIGGHAGDYGGEALDRTFQSFPTTKTYFPHFDLSPGSAQVKAHGKKVADALTTAVAHLDDLPGALSALSDLHAYKLRVDPVNFKLLSHCLLVTLACHHPTEFTPAVHASLDKFFTAVSTVLTSKYR"
    seq_gallus     = "MLTAEDKKLIQQAWEKAASHQEEFGAEALTRMFTTYPQTKTYFPHFDLSPGSDQVRGHGKKVLGALGNAVKNVDNLSQAMAELSNLHAYNLRVDPVNFKLLSQCIQVVLAVHMGKDYTPEVHAAFDKFLSAVSAVLAEKYR"
    seq_oncorhynchus = "XSLTAKDKSVVKAFWGKISGKADVVGAEALGRMLTAYPQTKTYFSHWADLSPGSGPVKKHGGIIMGAIGKAVGLMDDLVGGMSALSDLHAFKLRVDPGNFKILSHNILVTLAIHFPSDFTPEVHIAVDKFLAAVSAALADKYR"

    sequences = {
        "Chrysocyon brachyurus": seq_chrysocyon,
        "Gallus gallus": seq_gallus,
        "Oncorhynchus mykiss": seq_oncorhynchus
    }

    for name, seq in sequences.items():
        score, align1, align2 = smith_waterman(seq_hs, seq)
        print(f"--- Homo sapiens vs {name} ---")
        print(f"Maior pontuação: {score}")
        print(f"Subsequência alinhada (Homo sapiens): {align1}")
        print(f"Subsequência alinhada ({name}): {align2} \n")


if __name__ == "__main__":
    main()
