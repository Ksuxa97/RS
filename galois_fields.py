
class ExtendedGaloisField:
    def __init__(self, prime_order, n, primitive_poly):
        self.prime_order = prime_order
        self.primitive_poly = primitive_poly
        self.n = n
        self.log_bin, self.exp_bin = self.create_tables()

    def create_tables(self):

        binary = [0] * (self.n + 1)
        exp_binary = [0] * self.n * 2

        value = 1
        for i in range(self.n):
            exp_binary[i] = value
            exp_binary[i + self.n] = value
            binary[value] = i
            value = self.create_element(prev_element=value)

        return binary, exp_binary

    def bit_length(self, n):
        bits = 0
        while n >> bits: bits += 1
        return bits

    def create_element(self, prev_element):
        element = 0
        i = 0
        while (2 >> i) > 0:
            if 2 & (1 << i):
                element ^= prev_element << i
            i += 1

        dl1 = self.bit_length(element)
        dl2 = self.bit_length(self.primitive_poly)

        if dl1 < dl2:
            return element

        for i in range(dl1 - dl2, -1, -1):
            if element & (1 << i + dl2 - 1):
                element ^= self.primitive_poly << i
        return element

    def sum(self, x, y):
        return x ^ y

    def gf_mult(self, x, y):
        if x == 0 or y == 0:
            return 0

        return self.exp_bin[self.log_bin[x] + self.log_bin[y]]

    def gf_div(self, x, y):

        if x == 0 or y == 0:
            return 0

        power = (self.log_bin[x] - self.log_bin[y]) % self.n

        return self.exp_bin[power]

    def gf_pow(self, x, power):

        return self.exp_bin[(self.log_bin[x] * power) % self.n]

    def gf_inverse(self, x):
        return self.exp_bin[self.n - self.log_bin[x]]

