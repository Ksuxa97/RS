import random
import pandas
import numpy
from math import sqrt
from matplotlib import pyplot


from galois_fields import ExtendedGaloisField
from reed_solomon import RS_Extended
from irreducible_polynomials import irreducible_polynomials
from polinomials import Polynomials


def main():

    build_graph()

    return


###################################################
#        Functionality for graph creation         #
###################################################
def build_graph():

    create_statistics()
    read_csv("end_code.csv")


def create_statistics():
    m = 9
    b = 0
    prime_order = pow(2, m)
    n = prime_order - 1
    k = n // 2 + 1
    d = n - k + 1
    t = (d - 1) // 2
    primitive_poly = irreducible_polynomials.get(m)
    gf = ExtendedGaloisField(prime_order, n, primitive_poly)

    # generate root list
    roots = []
    for i in range(b, b + d - 1):
        roots.append(gf.exp_bin[i])

    poly = Polynomials(gf=gf)
    rs = RS_Extended(n=n, k=k, d=d, gf=gf, polynomial=poly, roots=roots)

    # calculate min poly
    generator = rs.generator_poly()
    result, check_g = rs.check_poly(generator)

    if not check_g:
        print("Something went wrong! Generator is not right.")

    print("G(x): ", generator)

    whole_list = []
    E_Ns = numpy.arange(0., 12., 0.5)
    sigmas = numpy.array([sqrt(1 / 10 ** (E_N / 10)) for E_N in E_Ns])

    counter = 0
    for sigma, E_N in zip(sigmas, E_Ns):
        for i in range(3000):
            info_sequence = create_info_sequence(n=n, k=k)
            c_x = rs.coder(generator=generator, info_sequence=info_sequence)
            y_x = canal_awgn(message=c_x, bits_for_number=m, sigma=sigma)

            wrong_word = 0
            decoded_message = rs.decoder(message=y_x, b=b, t=t)
            final_errors_count = compare_results(poly, c_x, decoded_message)

            if final_errors_count != 0:
                wrong_word = 1
            item = [
                sigma,
                E_N,
                final_errors_count,
                wrong_word
            ]
            whole_list.append(item)
            counter += 1

    np_array = numpy.array(whole_list)
    df = pandas.DataFrame({'Sigma': np_array[:, 0],
                           'E_N': np_array[:, 1],
                           'Final errors': np_array[:, 2],
                           'wrong_word': np_array[:, 3]
                           })

    name = 'end_code.csv'
    df.to_csv(name, index=False)

    return


def compare_results(polynomial, message, decoded_message):
    error_count = 0
    result = polynomial.bin_poly_sum(message, decoded_message)

    for r in result:
        if r != 0:
            error_count += 1

    return error_count


def read_csv(file_name):
    df = pandas.read_csv(file_name)
    result = {
        'E_b/N_0': [],
        'P': [],
    }
    for index in df['E_N'].unique():
        slice = df[df['E_N'] == index]
        result['E_b/N_0'].append(
            index
        )
        result['P'].append(
            slice['wrong_word'].sum() / 10
        )
    print(result['E_b/N_0'])
    print(result['P'])
    fig, ax = pyplot.subplots()
    ax.plot(result['E_b/N_0'], result['P'], '.', linestyle='-', linewidth=1, color='blue')
    ax.set_ylabel('p')
    ax.set_xlabel('E_b/N_0')
    pyplot.legend(loc='upper right')
    pyplot.yscale('log')
    pyplot.grid(True)


    pyplot.show()


###################################################
#  Testing encoding and decoding for m = 4....10  #
###################################################
def test():

    for m in range(4, 11):

        n = pow(2, m) - 1
        k = n // 2 + 1
        d = n - k + 1
        t = (d - 1) // 2
        b = 1

        print("M: ", m)

        error_count_t = 0
        error_count_t_1 = 0

        # generate field
        primitive_poly = irreducible_polynomials.get(m)
        gf = ExtendedGaloisField(prime_order=2, n=n, primitive_poly=primitive_poly)

        # generate root list
        roots = []
        for i in range(b, b + d - 1):
            roots.append(gf.exp_bin[i])

        poly = Polynomials(gf=gf)
        rs = RS_Extended(n=n, k=k, d=d, gf=gf, polynomial=poly, roots=roots)

        # calculate min poly
        generator = rs.generator_poly()
        result, check_g = rs.check_poly(generator)

        if not check_g:
            print("Something went wrong! Generator is not right.")

        print("G(x): ", generator)

        for j in range(10000):
            info_sequence = create_info_sequence(n=n, k=k)
            c_x = rs.coder(info_sequence=info_sequence, generator=generator)

            y_x1 = create_errors(c_x=c_x, t=t)
            y_x2 = create_errors(c_x=c_x, t=t + 1)

            res1 = calculate_errors(polynomials=poly, c_x=c_x, y_x=y_x1)
            res2 = calculate_errors(polynomials=poly, c_x=c_x, y_x=y_x2)

            decoded_msg1 = rs.decoder(message=y_x1, b=b, t=t)
            decoded_msg2 = rs.decoder(message=y_x2, b=b, t=t+1)

            error_count_t += compare_results(poly, c_x, decoded_msg1)
            error_count_t_1 += compare_results(poly, c_x, decoded_msg2)

        print("Error count for t: ", error_count_t)
        print("Error count for t+1: ", error_count_t_1)

        procent = error_count_t_1 / (((t + 1) * 10000) / 100)

        print("Probability: ", procent)


###################################################
#                  Simple runner                  #
###################################################
def run(m):
    b = 18
    n = pow(2, m) - 1
    k = 21
    d = n - k + 1
    t = (d - 1) // 2

    info_sequence = create_info_sequence(n=n, k=k)

    # generate field
    primitive_poly = irreducible_polynomials.get(m)
    gf = ExtendedGaloisField(2, m, primitive_poly)

    # generate root list
    roots = []
    for i in range(b, b + d - 1):
        roots.append(gf.log_bin[i%n])

    poly = Polynomials(gf=gf)
    rs = RS_Extended(n=n, k=k, d=d, gf=gf, polynomial=poly, roots=roots)

    c_x = rs.coder(info_sequence=info_sequence)

    y_x = canal_awgn(message=c_x, bits_for_number=m)

    # y = create_errors(c_x, t)

    res = calculate_errors(polynomials=poly, c_x=c_x, y_x=y_x)

    decoded_message = rs.decoder(message=y_x, b=b, t=t)

    return 0


def create_info_sequence(n, k):

    sequence = []

    for i in range(k):
        if i == k-1:
            symbol = random.randrange(1, n)
            sequence.append(symbol)
            break
        symbol = random.randrange(0, n)
        sequence.append(symbol)

    return sequence


def create_errors(c_x, t):

    y_x = c_x.copy()

    positions = []
    for i in range(len(y_x)):
        positions.append(i)

    for i in range(t):
        error_pos = random.choice(positions)
        value = y_x[error_pos]
        while value == y_x[error_pos]:
            value = random.randrange(0, len(y_x))
        y_x[error_pos] = value
        positions.remove(error_pos)

    return y_x


def canal_awgn(message, bits_for_number, sigma):

    y_bits = ''
    for i in range(len(message)):
        y_bits += "{0:b}".format(message[i]).zfill(bits_for_number)

    # convert bits to signal
    bits_list = list(y_bits)
    signal_array = [1. if i == '1' else -1. for i in bits_list]
    signal_array = numpy.array(signal_array)

    # add noise
    errors = numpy.random.normal(0, sigma, signal_array.shape[0])
    twisted_signal = signal_array + errors

    # convert signal to bits
    res = (twisted_signal >= 0).astype(int)
    bits = ''.join(str(i) for i in res)

    message_with_error = []
    for i in range(0, len(y_bits), bits_for_number):
        number = bits[i:i + bits_for_number]
        number = int(number, 2)
        message_with_error.append(number)

    return message_with_error


def calculate_errors(polynomials, c_x, y_x):

    poly = polynomials.bin_poly_sum(c_x, y_x)
    positions = []
    error_count = 0
    for i in range(len(poly)):
        if poly[i] != 0:
            error_count += 1
            positions.append(i)

    return poly, error_count, positions


if __name__ == '__main__':
    main()