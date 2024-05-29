import regex as re


def REpattern2regex(restriction_pattern):
	"""
	:param restriction_pattern: a restriction pattern
	:return: the regex formulation of the pattern
	"""
	if re.search("[^ACGTWN|\d]", restriction_pattern):
		raise ValueError("Error! invalid char is included")
	else:
		res = re.sub("N", "[ACGT]", restriction_pattern)
		res = re.sub("W", "[AT]", res)
		res = re.sub("(\d+)", r"{\1}", res)
		res = res.split("|")

	return res


def chop_by_restriction_enzymes(seq, restriction_pattern, additional_pattern=None):
	"""
	:param seq: a valid DNA sequence (i.e., over the alphabet ACGT, uppercase, and non-empty)
	:param restriction_pattern: a valid restriction enzyme's pattern, over the characters A, C, G, T, N (ambiguous for ACGT), and W (ambiguous for A or T)
	such that the restriction site is defined by a pipe "|", and a number that follows a character means that the character repeats that number of times
	:param additional_pattern: relevant for question 3.
	:return: a list of chops after applying the enzyme
	"""
	convert_pattern = REpattern2regex(restriction_pattern)
	prefix = convert_pattern[0]          # before the pipe
	suffix = convert_pattern[1]          # after the pipe
	if additional_pattern is None:
		pattern = f'(?<={prefix}){suffix}.*?{prefix}(?={suffix})'
		chops_lst = re.findall(pattern, seq)
	else:                  # additional_pattern is not None
		convert_additional_pattern = REpattern2regex(additional_pattern)
		additional_prefix = convert_additional_pattern[0]        # before the pipe
		additional_suffix = convert_additional_pattern[1]        # after the pipe
		pattern = f'(?<={prefix}){suffix}.*?{additional_prefix}(?={additional_suffix})'
		chops_lst = re.findall(pattern, seq)

	return chops_lst


def can_be_cloned(vector, insert, pattern1, pattern2):
	if pattern1 == pattern2:
		return False
	elif len(chop_by_restriction_enzymes(vector, pattern1, pattern2)) > 0 and len(chop_by_restriction_enzymes(insert, pattern1, pattern2)) > 0:
		return True
	elif len(chop_by_restriction_enzymes(vector, pattern2, pattern1)) > 0 and len(chop_by_restriction_enzymes(insert, pattern2, pattern1)) > 0:
		return True
	else:
		return False


def which_inserts_can_be_cloned (fasta_filepath, pattern1, pattern2):
	"""
	:param fasta_filepath: a path to a fasta file
	:param pattern1
	:param pattern2
	:return:
	"""
	open_file = open(fasta_filepath, 'r')
	read_lines = open_file.readlines()
	read_lines = [line.strip() for line in read_lines]

	vector_name = None
	vector_seq = ""
	insert_name = None
	insert_seq = ""
	inserts_can_be_cloned = []

	for line in read_lines:
		if line[0] == ">":             # new record
			if vector_name is None:
				vector_name = line[1:]      # name without >
			else:
				if insert_seq != "" and can_be_cloned(vector_seq, insert_seq, pattern1, pattern2):
					inserts_can_be_cloned.append(insert_name)
				insert_name = line[1:]         # update current insert name
				insert_seq = ""                # initialize insert sequence
		else:
			if insert_name is None:
				vector_seq += line
			else:
				insert_seq += line

	if can_be_cloned(vector_seq, insert_seq, pattern1, pattern2):     # check the last insert
		inserts_can_be_cloned.append(insert_name)

	return inserts_can_be_cloned


if __name__ == '__main__':
	# you can write whatever you want here
	pass
