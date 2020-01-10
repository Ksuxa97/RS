
def rs_farm(n, k, d, gf, polynomial, roots):

    if gf.prime_order == 2:
        rs = RS_Extended(n=n, k=k, d=d, gf=gf, polynomial= polynomial, roots=roots)

    else:
        rs = RS_Simple(gf=gf, polynomial=polynomial, roots=roots)

    return rs


class RS_Extended:

    def __init__(self, n, k, d, gf, polynomial, roots):
        self.gf = gf
        self.roots = roots
        self.polynomial = polynomial
        self.n = n
        self.k = k
        self.d = d
        self.p = n + 1

    def coder(self, info_sequence, generator):

        divider = [0] * (self.n - self.k)
        divider = divider + info_sequence

        quotient, reminder = self.polynomial.bin_poly_div(dividend=divider, divisor=generator)
        for i in range(len(reminder)):
            divider[i] = reminder[i]
        c = divider

        result, check_c = self.check_poly(c)

        if not check_c:
            print("The coding went wrong! C(x) is not correct.")

        print("C(x): ", c)

        return c

    def decoder(self, message, b, t):

        syndrome, all_zero_flag = self.check_poly(message)

        if all_zero_flag:
            return message

        for i in range(len(syndrome)-1, -1, -1):
            if syndrome[i] == 0:
                syndrome.pop()
            else:
                break

        # create x for Euclid
        x = [0] * (self.d - 1)
        x.append(1)

        error_value, error_locator = self.euclid(x=x, syndrome=syndrome, t=t)

        error_positions = self.find_error_positions(error_locator=error_locator)

        e_x = self.find_error_values(error_values=error_value, error_locator=error_locator,
                                     error_positions=error_positions)

        result = self.polynomial.bin_poly_sum(message, e_x)

        return result

    def generator_poly(self):

        g = [1]
        for root in self.roots:
            g = self.polynomial.bin_poly_mult(g, [root, 1])
        return g

    def check_poly(self, poly):

        result = []
        all_zero_flag = True

        for root in self.roots:
            value = self.calculate_polynomial(poly, root)
            if value != 0:
                all_zero_flag = False
            result.append(value)
        return result, all_zero_flag

    def euclid(self, x, syndrome, t):
        dividend = x.copy()
        divisor = syndrome.copy()

        quotient_list = []

        while True:
            quotient, reminder = self.polynomial.bin_poly_div(dividend, divisor)
            quotient_list.insert(0, quotient)
            if (len(reminder)-1) < t:
                break
            dividend = divisor
            divisor = reminder

        a = quotient_list[0]
        b = [1]
        for i in range(1, len(quotient_list)):
            temp = self.polynomial.bin_poly_mult(a, quotient_list[i])
            temp = self.polynomial.bin_poly_sum(temp, b)
            b = a
            a = temp

        return reminder, a

    def find_error_positions(self, error_locator):

        errors = []
        for i in range(self.n):
            root = self.gf.exp_bin[i]
            check = self.calculate_polynomial(error_locator, root)
            errors.append(check)

        error_positions = []
        for i in range(len(errors)):
            if errors[i] == 0:
                position = (self.n - i) % self.n
                error_positions.append(position)

        return error_positions

    def find_error_values(self, error_values, error_locator, error_positions):

        locator_derivative = self.calculate_derivative(poly=error_locator)

        # values = []
        e_x = [0] * self.n

        # res = self.alternative(error_values, error_locator, error_positions)

        for position in error_positions:
            root = self.gf.exp_bin[position]
            if position != 0:
                root = self.gf.exp_bin[self.n-position]

            locator_derivative_value = self.calculate_polynomial(locator_derivative, root)
            error_value = self.calculate_polynomial(error_values, root)
            value = self.gf.gf_div(error_value, locator_derivative_value)
            e_x[position] = value

        return e_x

    def alternative_forney(self, error_values, error_positions, b):

        values = [0] * self.n

        for error_position in error_positions:
            x_j = error_positions.copy()
            x_j.remove(error_position)
            x_l = self.gf.log_bin[(self.n-error_position)%self.n]
            x_l_b = self.gf.gf_pow(self.gf.log_bin[error_position], -b)

            devisor = 1
            for i in range(len(x_j)):
                temp = self.gf.log_bin[x_j[i]]
                temp = 1 ^ self.gf.gf_mult(temp, x_l)
                devisor = self.gf.gf_mult(devisor, temp)

            error_value = self.calculate_polynomial(error_values, x_l)
            dividend = self.gf.gf_mult(error_value, x_l_b)

            result = self.gf.gf_div(dividend, devisor)

            values[error_position] = result

        return values

    def calculate_derivative(self, poly):

        derivative = []

        for i in range(1, len(poly)):
            if i % 2 == 0 and i == len(poly)-1:
                break
            elif i % 2 == 0:
                derivative.append(0)
            else:
                derivative.append(poly[i])

        return derivative

    def calculate_polynomial(self, poly, root):

        result = 0

        for i in range(len(poly)):
            temp_res = self.gf.gf_pow(root, i)
            temp_res = self.gf.gf_mult(temp_res, poly[i])
            result = self.gf.sum(result, temp_res)

        return result


class RS_Simple:
    def __init__(self, gf, polynomial, roots):
        self.gf = gf
        self.roots = roots
        self.polynomial = polynomial

    def generator_poly(self):

        g = [1]
        for root in self.roots:
            g = self.polynomial.poly_mult(g, [1, self.gf.prime_order - root])
        return g

    def check_poly(self, poly):

        result = []

        temp_poly = poly.copy()
        temp_poly.reverse()

        for root in self.roots:
            temp_res = 0
            for i in range(len(temp_poly)):
                element = self.gf.gf_pow(root, i)
                element = self.gf.gf_mult(element, temp_poly[i])
                temp_res = self.gf.sum(temp_res, element)
            result.insert(0, temp_res)
        return result

    def euclid(self, x, syndrome, t):
        dividend = x.copy()
        divisor = syndrome.copy()
        quotient = 0
        reminder = 0

        while True:
            quotient, reminder = self.polynomial.poly_div(dividend, divisor)
            if len(reminder) <= t+1:
                break
            dividend = divisor
            divisor = reminder

        return quotient, reminder, divisor