from datetime import datetime
from itertools import product
from numpy import array, longlong, min, max, where
from warnings import filterwarnings

filterwarnings("ignore", category=RuntimeWarning)


def encode(binary_message, strand_index, pattern, mapping, bio_filter,
           salt_number=46, previous_number=8, low_order_number=10):
    """
    Encode the binary message.

    :param binary_message: binary message.
    :type binary_message: numpy.ndarray
    :param strand_index: index of current strand.
    :type strand_index: int
    :param pattern: pattern in HEDGES paper (Table 3).
    :type pattern: list
    :param mapping: mapping between 0/1/2/3 and A/C/G/T.
    :type mapping: list
    :param bio_filter: established local constraint set.
    :type bio_filter: LocalBioFilter
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
    """
    dna_string, available_nucleotides, bit_location, pattern_flag = "", mapping, 0, 0
    salt_index = strand_index % (2 ** salt_number)  # s bits of salt (strand ID).
    while bit_location < len(binary_message):
        bit_index = bit_location % (2 ** low_order_number)  # low-order q bits of the bit position index i.

        if bit_location - previous_number >= 0:
            previous_info = binary_message[bit_location - previous_number: bit_location]
            previous_value = bit_to_number(previous_info)
        else:
            previous_value = 0

        hash_value = hash_function(bit_index | previous_value | salt_index)
        if len(available_nucleotides) > 0:
            bit_number = pattern[pattern_flag]
            message_bit = binary_message[bit_location: bit_location + bit_number]
            bit_value = bit_to_number(message_bit) if len(message_bit) > 0 else 0
            nucleotide = available_nucleotides[(hash_value + bit_value) % len(available_nucleotides)]
            bit_location += bit_number
            pattern_flag = (pattern_flag + 1) % len(pattern)
        else:
            raise ValueError("DNA string (index = " + str(strand_index) + ") " +
                             "cannot be encoded because of the established constraints!")

        dna_string += nucleotide

        available_nucleotides = []
        for potential_nucleotide in mapping:
            if bio_filter.valid(dna_string + potential_nucleotide, only_last=True):
                available_nucleotides.append(potential_nucleotide)

    return dna_string


def decode(dna_string, strand_index, bit_length, pattern, mapping, bio_filter,
           salt_number=46, previous_number=8, low_order_number=10, heap_limitation=1e6,
           initial_score=0.0, correct_penalty=-0.035, insert_penalty=1.0, delete_penalty=1.0, mutate_penalty=1.0):
    """
    Decode or repair the DNA string based on A star search.

    :param dna_string: encoded DNA string.
    :type dna_string: str
    :param strand_index: index of current strand.
    :type strand_index: int
    :param bit_length: length of binary message.
    :type bit_length: int
    :param pattern: pattern in HEDGES paper (Table 3).
    :type pattern: list
    :param mapping: mapping between 0/1/2/3 and A/C/G/T.
    :type mapping: list
    :param bio_filter: established local constraint set.
    :type bio_filter: LocalBioFilter
    :param salt_number: number of salt bits.
    :type salt_number: int
    :param previous_number: number of previous bits.
    :type previous_number: int
    :param low_order_number: number of low-order bits.
    :type low_order_number: int
    :param heap_limitation: limitation of the size of heap.
    :type heap_limitation: int
    :param initial_score: initial score starting the A* search.
    :type initial_score: float
    :param correct_penalty: penalty when correction.
    :type correct_penalty: float
    :param insert_penalty: penalty when insertion.
    :type insert_penalty: float
    :param delete_penalty: penalty when deletion.
    :type delete_penalty: float
    :param mutate_penalty: penalty when mutation.
    :type mutate_penalty: float

    :return: decoded or repaired information.
    :rtype: (list, int, int)

    .. note::
        The returned information includes three parts.
        The first element is a pair list containing the won candidates,
        the pair of which is a binary message and its corresponding DNA string.
        The second element is the final size of heap.
        The third element is the reading process of DNA string (after search).
    """
    # s bits of salt (strand index).
    salt_value = strand_index % (2 ** salt_number)

    class HypothesisNode:

        def __init__(self, pattern_flag, message, string):
            self.pattern_flag, self.message, self.string = pattern_flag, message, string

        def next(self, nucleotide_index, current_score):
            follow_vertices, follow_scores, follow_indices = [], [], []

            # collect the available nucleotides in this location.
            available_nucleotides = []
            for potential_nucleotide in mapping:
                if bio_filter.valid(self.string + potential_nucleotide, only_last=True):
                    available_nucleotides.append(potential_nucleotide)

            if len(available_nucleotides) == 0:  # this path is blocked, stop running.
                return [], [], [], []

            # low-order q bits of the bit position index i.
            bit_index = len(self.message) % (2 ** low_order_number)

            if len(self.message) - previous_number >= 0:  # p previous concatenated bits.
                previous_info = self.message[len(self.message) - previous_number:]
                previous_value = bit_to_number(previous_info)
            else:
                previous_value = 0

            for message_bit in product([0, 1], repeat=pattern[self.pattern_flag]):
                hash_value = hash_function(bit_index | previous_value | salt_value)
                bit_value = bit_to_number(list(message_bit)) if len(message_bit) > 0 else 0
                nucleotide = available_nucleotides[(hash_value + bit_value) % len(available_nucleotides)]
                message, string = self.message + list(message_bit), self.string + nucleotide

                if nucleotide == dna_string[nucleotide_index]:  # assume that current nucleotide is correct.
                    follow_vertices.append(HypothesisNode((self.pattern_flag + 1) % len(pattern), message, string))
                    follow_scores.append(current_score + correct_penalty)
                    follow_indices.append(nucleotide_index + 1)
                else:
                    # assume that current nucleotide is mutated.
                    follow_vertices.append(HypothesisNode((self.pattern_flag + 1) % len(pattern), message, string))
                    follow_scores.append(current_score + mutate_penalty)
                    follow_indices.append(nucleotide_index + 1)

                    # assume that current nucleotide is inserted, the (i + 1)-th nucleotide is i-th nucleotide.
                    if nucleotide_index + 1 < len(dna_string) and nucleotide == dna_string[nucleotide_index + 1]:
                        follow_vertices.append(HypothesisNode((self.pattern_flag + 1) % len(pattern), message, string))
                        follow_scores.append(current_score + insert_penalty)
                        follow_indices.append(nucleotide_index + 2)

                    # assume that current nucleotide is deleted.
                    follow_vertices.append(HypothesisNode((self.pattern_flag + 1) % len(pattern), message, string))
                    follow_scores.append(current_score + delete_penalty)
                    follow_indices.append(nucleotide_index)

            return follow_vertices, follow_scores, follow_indices, [len(v.message) for v in follow_vertices]

    monitor, terminal_indices, heap_size = Monitor(), None, 1
    heap = {"v": [HypothesisNode(0, [], "")], "s": [initial_score], "i": [0], "l": [0]}  # priority heap

    while True:  # repair by A star search (score priority).
        chuck_indices, chuck_score = where(array(heap["s"]) == min(heap["s"]))[0], min(heap["s"])

        for chuck_index in chuck_indices:
            # noinspection PyTypeChecker
            heap["s"][chuck_index] = len(dna_string) + 1  # set the chuck vertex to inaccessible.
            heap_size -= 1

            if heap["i"][chuck_index] >= len(dna_string):
                continue

            follow_info = heap["v"][chuck_index].next(heap["i"][chuck_index], chuck_score)
            heap_size += len(follow_info[0])

            heap["v"], heap["l"] = heap["v"] + follow_info[0], heap["l"] + follow_info[3]
            heap["s"], heap["i"] = heap["s"] + follow_info[1], heap["i"] + follow_info[2]

            current_process = heap["i"][chuck_index] + 1
            monitor.output(current_process, len(dna_string), extra={"size": len(heap["v"]), "score": chuck_score})

            # the first chain of hypotheses to decode the required bytes of message wins.
            if bit_length == max(heap["l"]) or len(heap["v"]) >= heap_limitation:
                if current_process < len(dna_string):
                    print()  # not finish.

                candidates = []
                for terminal_index in where(array(heap["l"]) == bit_length)[0]:
                    candidates.append((heap["v"][terminal_index].message, heap["v"][terminal_index].string))

                return candidates, heap_size, current_process


def hash_function(source_value):
    """
    Obtain the target value from the source value based on the well-accepted hash function.

    :param source_value: source bit value.
    :type source_value: int

    :return: target value after the hash function.
    :rtype: int
    """
    target_value = array(source_value, dtype=longlong) * array(3935559000370003845, dtype=longlong)
    target_value += array(2691343689449507681, dtype=longlong)
    target_value ^= target_value >> 21
    target_value ^= target_value << 37
    target_value ^= target_value >> 4
    target_value *= array(4768777513237032717, dtype=longlong)
    target_value ^= target_value << 20
    target_value ^= target_value >> 41
    target_value ^= target_value << 5

    return target_value


def bit_to_number(bit_array):
    """
    Transform a bit array to the equivalent decimal number.

    :param bit_array: bit array.
    :type bit_array: list

    :return: equivalent decimal number (may huge) of the inputted bit array.
    :rtype: str or int
    """
    decimal_number = 0

    for a_bit in bit_array:
        decimal_number = decimal_number * 2 + a_bit

    return decimal_number


class LocalBioFilter(object):

    def __init__(self, observed_length, max_homopolymer_runs=None, gc_range=None, undesired_motifs=None):
        """
        Initialize the screen of local biochemical constraints.

        :param observed_length: length of the DNA string observed in the window.
        :type observed_length: int

        :param max_homopolymer_runs: maximum homopolymer runs.
        :type max_homopolymer_runs: int

        :param gc_range: range of GC content.
        :type gc_range: list

        :param undesired_motifs: undesired DNA motifs.
        :type undesired_motifs: list
        """
        self.screen_name = "Local"
        if max_homopolymer_runs is not None:
            if observed_length < max_homopolymer_runs:
                raise ValueError("The parameter \"observed_length\" must "
                                 + "longer than the parameter \"max_homopolymer_runs\"!")
        if undesired_motifs is not None:
            for index, undesired_motif in enumerate(undesired_motifs):
                if len(undesired_motif) > observed_length:
                    raise ValueError("The parameter \"observed_length\" must "
                                     + "longer than the length of any motif in the parameter \"undesired_motifs\"!")

        self._observed_length = observed_length
        self._max_homopolymer_runs = max_homopolymer_runs
        self._gc_range = gc_range
        self._undesired_motifs = undesired_motifs

    def valid(self, dna_string, only_last=True):
        """
        Judge whether the DNA string meets the local biochemical constraints.

        :param dna_string: DNA string to be judged.
        :type dna_string: str

        :param only_last: only check the DNA string of the last observed window.
        :type only_last: bool

        :return: judgement.
        :rtype: bool

        .. note::
            "only_last" parameter is used to save time.
            For most tree-based coding algorithms,
            it is not necessary to detect the sub DNA strings observed in each window from scratch every time.
        """
        if only_last:
            observed_dna_string = dna_string[-self._observed_length:]
        else:
            observed_dna_string = dna_string

        for nucleotide in observed_dna_string:
            if nucleotide not in "ACGT":
                return False

        if self._max_homopolymer_runs is not None:
            for nucleotide in "ACGT":
                if nucleotide * (1 + self._max_homopolymer_runs) in observed_dna_string:
                    return False

        if self._undesired_motifs is not None:
            for special in self._undesired_motifs:
                if special in observed_dna_string:
                    return False
                reverse_complement = special.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
                reverse_complement = reverse_complement[::-1].upper()
                if reverse_complement in observed_dna_string:
                    return False

        if self._gc_range is not None:
            if len(observed_dna_string) >= self._observed_length:
                for index in range(len(observed_dna_string) - self._observed_length + 1):
                    sub_dna_string = observed_dna_string[index: index + self._observed_length]
                    gc_count = sub_dna_string.count("C") + sub_dna_string.count("G")
                    if gc_count > self._gc_range[1] * self._observed_length:
                        return False
                    if gc_count < self._gc_range[0] * self._observed_length:
                        return False
            else:
                gc_count = observed_dna_string.count("C") + observed_dna_string.count("G")
                if gc_count > self._gc_range[1] * self._observed_length:
                    return False
                at_count = observed_dna_string.count("A") + observed_dna_string.count("T")
                if at_count > (1 - self._gc_range[0]) * self._observed_length:
                    return False

        return True

    def __str__(self):
        info = self.screen_name + "\n"
        info += "maximum homopolymer runs : " + str(self._max_homopolymer_runs) + "\n"
        info += "local GC content range   : " + str(self._gc_range[0]) + "<= GC <=" + str(self._gc_range[1]) + "\n"
        info += "undesired DNA motifs     : " + str(self._undesired_motifs).replace("\"", "") + "\n"
        return info


class Monitor(object):

    def __init__(self):
        """
        Initialize the monitor to identify the task progress.
        """
        self.last_time = None

    def output(self, current_state, total_state, extra=None):
        """
        Output the current state of process.

        :param current_state: current state of process.
        :type current_state: int

        :param total_state: total state of process.
        :type total_state: int

        :param extra: extra vision information if required.
        :type extra: dict
        """
        if self.last_time is None:
            self.last_time = datetime.now()

        if current_state == 0:
            return

        position = int(current_state / total_state * 100)

        string = "|"

        for index in range(0, 100, 5):
            if position >= index:
                string += "█"
            else:
                string += " "

        string += "|"

        pass_time = (datetime.now() - self.last_time).total_seconds()
        wait_time = int(pass_time * (total_state - current_state) / current_state)

        string += " " * (3 - len(str(position))) + str(position) + "% ("

        string += " " * (len(str(total_state)) - len(str(current_state))) + str(current_state) + "/" + str(total_state)

        if current_state < total_state:
            minute, second = divmod(wait_time, 60)
            hour, minute = divmod(minute, 60)
            string += ") wait " + "%04d:%02d:%02d" % (hour, minute, second)
        else:
            minute, second = divmod(pass_time, 60)
            hour, minute = divmod(minute, 60)
            string += ") used " + "%04d:%02d:%02d" % (hour, minute, second)

        if extra is not None:
            string += " " + str(extra).replace("\'", "").replace("{", "(").replace("}", ")") + "."
        else:
            string += "."

        print("\r" + string, end="", flush=True)

        if current_state >= total_state:
            self.last_time = None
            print()
