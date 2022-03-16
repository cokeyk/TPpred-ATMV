def check_rewrite_data(data, data_path):
    with open(data_path, 'w') as output_data:
        first_is_name = False
        need_line = False
        seq_num = 0
        for x in data:
            x = x.strip()
            if not x:
                continue
            if x[0] == '>':
                x = adapt_PEPformat(x)      # for PEPred-Suite software
                seq_num += 1
                if need_line:
                    return False, 'Sorry, sequence ' + str(seq_num - 1) + ' is empty, please check it!'
                need_line = True
                first_is_name = True
            else:
                if not first_is_name:
                    seq_num += 1
                    return False, 'Sorry, the name of sequence ' + str(seq_num) + ' is wrong, please check it!'
                if not is_protein_seq(x):
                    return False, 'Sorry, sequence ' + str(seq_num) + ' contains wrong alphabet, please check it!'
                if not need_line:
                    seq_num += 1
                    return False, 'Sorry, the name of sequence ' + str(seq_num) + ' is wrong, please check it!'
                need_line = False
            output_data.write(x + '\n')
        if need_line:
            return False, 'Sorry, sequence ' + str(seq_num) + ' is empty, please check it!'
        return True, ''

def adapt_PEPformat(x):
    """
    PEPred-Suite jar input format
    :param x: id that starts with '>'
    :return: id that starts with '>' and ends with '|0'
    """
    if '|' not in x:
        x = x + '|0'
    return x

def is_protein_seq(seq):
    nchar = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    for x in seq:
        if x not in nchar:
            return False
    return True