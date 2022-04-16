from itertools import product
from numpy import array, longlong, floor, log2, max, where
from warnings import filterwarnings
from dsw import bit_to_number, Monitor

filterwarnings("ignore", category=RuntimeWarning)


def hash_function(source_value):
    """
    Obtain the target value from the source value based on the well-accepted hash function.

    :param source_value: source bit value.
    :type source_value: int

    :return: target value after the hash function.
    :rtype: int
    """
    target_value = source_value * array(3935559000370003845, dtype=longlong)
    target_value += array(2691343689449507681, dtype=longlong)
    target_value ^= target_value >> 21
    target_value ^= target_value << 37
    target_value ^= target_value >> 4
    target_value *= array(4768777513237032717, dtype=longlong)
    target_value ^= target_value << 20
    target_value ^= target_value >> 41
    target_value ^= target_value << 5

    return target_value


def encode(binary_message, strand_index, mapping, bio_filter,
           salt_number=46, previous_number=8, low_order_number=10):
    """
    Encode the binary message.

    :param binary_message: binary message.
    :type binary_message: numpy.ndarray
    :param strand_index: index of current strand.
    :type strand_index: int
    :param mapping: mapping between 0/1/2/3 and A/C/G/T.
    :type mapping: list
    :param bio_filter: established local constraint set.
    :type bio_filter: dsw.LocalBioFilter
    :param salt_number: number of salt bits.
    :type salt_number: int
    :param previous_number: number of previous bits.
    :type previous_number: int
    :param low_order_number: number of low-order bits.
    :type low_order_number: int

    :return: encoded DNA string.
    :rtype: str

    .. note::
        The parameter "mapping" is the order of A/C/G/T.
        For example, if "mapping" is [C, G, T, A], the mapping between digit and nucleotide is
        0-C, 1-G, 2-T, and 3-T.

        This is an extended version that can accept any local biochemical constraints.
        The treatment when the number of available nucleotides is 3 is not clearly described in Press' work.
        Here, when it equals 3, we only take the first two available nucleotides.
    """
    dna_string, available_nucleotides, bit_location = "", mapping, 0
    salt_index = strand_index % (2 ** salt_number)  # s bits of salt (strand ID).
    while bit_location < len(binary_message):
        bit_index = bit_location % (2 ** low_order_number)  # low-order q bits of the bit position index i.

        if bit_location - previous_number >= 0:
            previous_info = binary_message[bit_location - previous_number: bit_location]
            previous_value = bit_to_number(previous_info, is_string=False)
        else:
            previous_value = 0

        if len(available_nucleotides) == 1:
            nucleotide = available_nucleotides[0]

        elif len(available_nucleotides) == 2 or len(available_nucleotides) == 3:
            hash_value = hash_function(bit_index | previous_value | salt_index) % 2
            bit_value = binary_message[bit_location]
            nucleotide = available_nucleotides[(hash_value + bit_value) % 2]
            bit_location += 1

        else:
            hash_value = hash_function(bit_index | previous_value | salt_index) % 4
            if bit_location + 2 <= len(binary_message):
                bit_value = binary_message[bit_location] * 2 + binary_message[bit_location + 1]
            else:
                bit_value = binary_message[bit_location]
            nucleotide = available_nucleotides[(hash_value + bit_value) % 4]
            bit_location += 2

        dna_string += nucleotide

        available_nucleotides = []
        for potential_nucleotide in mapping:
            if bio_filter.valid(dna_string + potential_nucleotide, only_last=True):
                available_nucleotides.append(potential_nucleotide)

        if len(available_nucleotides) == 0:
            raise ValueError("DNA string (index = " + str(strand_index) + ") " +
                             "cannot be encoded because of the established constraints!")

    return dna_string


def decode(dna_string, strand_index, bit_length, mapping, bio_filter,
           salt_number=46, previous_number=8, low_order_number=10):
    """
    Decode the DNA string.

    :param dna_string: encoded DNA string.
    :type dna_string: str
    :param strand_index: index of current strand.
    :type strand_index: int
    :param bit_length: length of binary message.
    :type bit_length: int
    :param mapping: mapping between 0/1/2/3 and A/C/G/T.
    :type mapping: list
    :param bio_filter: established local constraint set.
    :type bio_filter: dsw.LocalBioFilter
    :param salt_number: number of salt bits.
    :type salt_number: int
    :param previous_number: number of previous bits.
    :type previous_number: int
    :param low_order_number: number of low-order bits.
    :type low_order_number: int

    :return: decoded binary message.
    :rtype: numpy.ndarray

    .. note::
        The parameter "mapping" is the order of A/C/G/T.
        For example, if "mapping" is [C, G, T, A], the mapping between digit and nucleotide is
        0-C, 1-G, 2-T, and 3-T.

        This is an extended version that can accept any local biochemical constraints.
        The treatment when the number of available nucleotides is 3 is not clearly described in Press' work.
        Here, when it equals 3, we only take the first two available nucleotides.
    """
    binary_message, available_nucleotides = [], mapping
    salt_index = strand_index % (2 ** salt_number)  # s bits of salt (strand ID).

    for nucleotide_index, nucleotide in enumerate(dna_string):
        bit_index = len(binary_message) % (2 ** low_order_number)  # low-order q bits of the bit position index i.
        if len(binary_message) - previous_number >= 0:
            previous_info = binary_message[len(binary_message) - previous_number:]
            previous_value = bit_to_number(previous_info, is_string=False)
        else:
            previous_value = 0

        if len(available_nucleotides) == 1:
            pass

        elif len(available_nucleotides) == 2 or len(available_nucleotides) == 3:
            hash_value = hash_function(bit_index | previous_value | salt_index) % 2
            if available_nucleotides[(hash_value + 0) % 2] == nucleotide:
                binary_message.append(0)
            else:
                binary_message.append(1)
        else:
            hash_value = hash_function(bit_index | previous_value | salt_index) % 4
            for bit_value in range(4):
                if available_nucleotides[(hash_value + bit_value) % 4] == nucleotide:
                    if len(binary_message) + 2 > bit_length:
                        binary_message.append(bit_value % 2)
                    else:
                        binary_message += [int(bit_value / 2), bit_value % 2]
                    break

        available_nucleotides = []
        for potential_nucleotide in mapping:
            if bio_filter.valid(dna_string[:nucleotide_index + 1] + potential_nucleotide, only_last=True):
                available_nucleotides.append(potential_nucleotide)

        if len(available_nucleotides) == 0:
            raise ValueError("DNA string (index = " + str(strand_index) + ") contains error(s)!")

    binary_message = binary_message[:bit_length]

    return array(binary_message)


def repair(dna_string, strand_index, initial_score, bit_length, mapping, bio_filter,
           salt_number=46, previous_number=8, low_order_number=10, heap_limitation=1e6,
           correct_penalty=-0.035, insert_penalty=1.0, delete_penalty=1.0, mutate_penalty=1.0):
    """
    Repair the wrong DNA string based on A* search.

    :param dna_string: encoded DNA string.
    :type dna_string: str
    :param strand_index: index of current strand.
    :type strand_index: int
    :param initial_score: initial score starting the A* search.
    :type initial_score: float
    :param bit_length: length of binary message.
    :type bit_length: int
    :param mapping: mapping between 0/1/2/3 and A/C/G/T.
    :type mapping: list
    :param bio_filter: established local constraint set.
    :type bio_filter: dsw.LocalBioFilter
    :param salt_number: number of salt bits.
    :type salt_number: int
    :param previous_number: number of previous bits.
    :type previous_number: int
    :param low_order_number: number of low-order bits.
    :type low_order_number: int
    :param heap_limitation: limitation of the size of heap.
    :type heap_limitation: int
    :param correct_penalty: penalty when correction.
    :type correct_penalty: float
    :param insert_penalty: penalty when insertion.
    :type insert_penalty: float
    :param delete_penalty: penalty when deletion.
    :type delete_penalty: float
    :param mutate_penalty: penalty when mutation.
    :type mutate_penalty: float

    :return: repaired DNA string set (if not one solution), heap limitation
    :rtype: (list, int)

    .. note::
        It is not the original design.

        Here, 6 child hypotheses are expanded to more child hypotheses
        to suppose the unknown number of bits in the current message bit.
    """

    class InfoVertex:

        def __init__(self, previous_value, bit_index, string):
            self.previous_value, self.bit_index, self.string = previous_value, bit_index, string

        # 6 child hypotheses are expanded to more child hypotheses (at most under the established constraints).
        def next(self, salt_value, nucleotide_index, current_score):
            follow_vertices, follow_scores, follow_indices = [], [], []

            if nucleotide_index > len(dna_string) - 1:
                return [], [], [], []

            # collect the available nucleotides in this location.
            available_nucleotides = []
            for potential_nucleotide in mapping:
                if bio_filter.valid(self.string + potential_nucleotide, only_last=True):
                    available_nucleotides.append(potential_nucleotide)

            if len(available_nucleotides) == 0:  # this path is blocked, stop running.
                return [], [], [], []

            bit_index = self.bit_index % (2 ** low_order_number)  # low-order q bits of the bit position index i.

            # assume that current nucleotide is correct or is mutated.
            if dna_string[nucleotide_index] in available_nucleotides:
                nucleotide = dna_string[nucleotide_index]
                for message_bit in product([0, 1], repeat=int(floor(log2(len(available_nucleotides))))):
                    if len(message_bit) == 1:
                        hash_value = hash_function(bit_index | self.previous_value | salt_value) % 2
                        previous_value = (self.previous_value * 2 + message_bit[0]) % (2 ** previous_number)
                        bit_value = message_bit[0]
                        if available_nucleotides[(hash_value + bit_value) % 2] == nucleotide:
                            vertex = InfoVertex(previous_value, self.bit_index + 1, self.string + nucleotide)
                            score = current_score + correct_penalty
                        else:
                            vertex = InfoVertex(previous_value, self.bit_index + 1, self.string + nucleotide)
                            score = current_score + mutate_penalty
                    elif len(message_bit) == 2:
                        hash_value = hash_function(bit_index | self.previous_value | salt_value) % 4
                        bit_value = message_bit[0] * 2 + message_bit[1]
                        previous_value = (self.previous_value * 4 + bit_value) % (2 ** previous_number)
                        if available_nucleotides[(hash_value + bit_value) % 4] == nucleotide:
                            vertex = InfoVertex(previous_value, self.bit_index + 2, self.string + nucleotide)
                            score = current_score + correct_penalty
                        else:
                            vertex = InfoVertex(previous_value, self.bit_index + 2, self.string + nucleotide)
                            score = current_score + mutate_penalty
                    else:  # len(message_bit) == 0.
                        vertex = InfoVertex(self.previous_value, self.bit_index, self.string + nucleotide)
                        score = current_score + correct_penalty
                    follow_vertices.append(vertex)
                    follow_scores.append(score)
                    follow_indices.append(nucleotide_index + 1)

            # assume that the current nucleotide is inserted
            # so that the right current nucleotide (i) is its latter nucleotide (i + 1).
            if nucleotide_index + 1 < len(dna_string):
                if dna_string[nucleotide_index + 1] in available_nucleotides:  # move i to i + 1.
                    nucleotide = dna_string[nucleotide_index + 1]
                    for message_bit in product([0, 1], repeat=int(floor(log2(len(available_nucleotides))))):
                        vertex = None
                        if len(message_bit) == 1:
                            hash_value = hash_function(bit_index | self.previous_value | salt_value) % 2
                            bit_value = message_bit[0]
                            previous_value = (self.previous_value * 2 + message_bit[0]) % (2 ** previous_number)
                            if available_nucleotides[(hash_value + bit_value) % 2] == nucleotide:
                                vertex = InfoVertex(previous_value, self.bit_index + 1, self.string + nucleotide)
                        elif len(message_bit) == 2:
                            hash_value = hash_function(bit_index | self.previous_value | salt_value) % 4
                            bit_value = message_bit[0] * 2 + message_bit[1]
                            previous_value = (self.previous_value * 4 + bit_value) % (2 ** previous_number)
                            if available_nucleotides[(hash_value + bit_value) % 4] == nucleotide:
                                vertex = InfoVertex(previous_value, self.bit_index + 2, self.string + nucleotide)
                        else:  # len(message_bit) == 0.
                            vertex = InfoVertex(self.previous_value, self.bit_index, self.string + nucleotide)
                        if vertex is not None:
                            follow_vertices.append(vertex)
                            follow_scores.append(current_score + insert_penalty)
                            follow_indices.append(nucleotide_index + 2)

            # assume that current nucleotide is deleted
            # so that the right current nucleotide (i) is an available nucleotide.
            for message_bit in product([0, 1], repeat=int(floor(log2(len(available_nucleotides))))):
                if len(message_bit) == 1:
                    hash_value = hash_function(bit_index | self.previous_value | salt_value) % 2
                    bit_value = message_bit[0]
                    previous_value = (self.previous_value * 2 + bit_value) % (2 ** previous_number)
                    nucleotide = available_nucleotides[(hash_value + bit_value) % 2]
                    vertex = InfoVertex(previous_value, self.bit_index + 1, self.string + nucleotide)
                    follow_vertices.append(vertex)
                    follow_scores.append(current_score + delete_penalty)
                    follow_indices.append(nucleotide_index)
                elif len(message_bit) == 2:
                    hash_value = hash_function(bit_index | self.previous_value | salt_value) % 4
                    bit_value = message_bit[0] * 2 + message_bit[1]
                    previous_value = (self.previous_value * 4 + bit_value) % (2 ** previous_number)
                    nucleotide = available_nucleotides[(hash_value + message_bit[0]) % 4]
                    vertex = InfoVertex(previous_value, self.bit_index + 2, self.string + nucleotide)
                    follow_vertices.append(vertex)
                    follow_scores.append(current_score + delete_penalty)
                    follow_indices.append(nucleotide_index)
                else:  # len(message_bit) == 0 (available nucleotide choices are NOT be limited).
                    for nucleotide in available_nucleotides:
                        vertex = InfoVertex(self.previous_value, self.bit_index, self.string + nucleotide)
                        follow_vertices.append(vertex)
                        follow_scores.append(current_score + delete_penalty)
                        follow_indices.append(nucleotide_index)

            return follow_vertices, follow_scores, follow_indices, [v.bit_index for v in follow_vertices]

    monitor, terminal_indices = Monitor(), None
    heap = {"v": [InfoVertex(0, 0, "")], "s": [initial_score], "i": [0], "l": [0]}  # priority heap

    # repair by A star search (score priority).
    while True:
        used_indices = where(array(heap["s"]) == min(heap["s"]))[0]
        used_score = heap["s"][used_indices[0]]
        finished = False
        for chuck in used_indices:
            # set the chuck vertex to inaccessible.
            used_vertex, base_index, heap["s"][chuck] = heap["v"][int(chuck)], heap["i"][int(chuck)], len(dna_string)
            follow_info = used_vertex.next(strand_index % (2 ** salt_number), base_index, used_score)

            monitor.output(base_index, len(dna_string), extra={"size": len(heap["v"]), "score": used_score})

            # ignore the limitation of heap size.
            # the first chain of hypotheses to decode the required bytes of message wins.
            if bit_length == max(heap["l"]) or len(heap["v"]) > heap_limitation:
                lengths = array(heap["l"])
                terminal_indices = where(lengths == bit_length)[0]
                finished = True
                break

            heap["v"], heap["l"] = heap["v"] + follow_info[0], heap["l"] + follow_info[3]
            heap["s"], heap["i"] = heap["s"] + follow_info[1], heap["i"] + follow_info[2]

        if finished:
            results = []
            for terminal_index in terminal_indices:
                # DNA string with a redundancy nucleotide at last.
                results.append(heap["v"][terminal_index].string)
            return list(set(results)), len(heap["v"])
