class Polynomials:
    def __init__(self, gf):
        self.gf = gf

    def bin_poly_sum(self, x, y):
        x_temp = x.copy()
        y_temp = y.copy()

        if len(x_temp) > len(y_temp):
            for i in range(len(x_temp) - len(y_temp)):
                y_temp.append(0)
        elif len(y_temp) > len(x_temp):
            for i in range(len(y_temp) - len(x_temp)):
                x_temp.append(0)

        sum = []
        for i in range(len(x_temp)):
            temp_sum = self.gf.sum(x_temp[i], y_temp[i])
            sum.append(temp_sum)

        return sum

    def bin_poly_mult(self, p, q):

        r = [0] * (len(p) + len(q) - 1)

        for j in range(len(q)-1, -1, -1):
            for i in range(len(p)-1, -1, -1):
                r[i + j] ^= self.gf.gf_mult(p[i], q[j])
        return r

    def bin_poly_div(self, dividend, divisor):

        msg_out = list(dividend)
        quotient = [0] * (len(dividend) - len(divisor) + 1)

        while len(msg_out) >= len(divisor):

            i = len(msg_out) - (len(divisor) - 1) - 1

            coef = self.gf.gf_div(msg_out[len(msg_out) - 1], divisor[len(divisor) - 1])
            quotient[i] = coef
            if coef != 0:

                for j in range(len(divisor)-1, -1, -1):
                    if divisor[j] != 0:
                        msg_out[i + j] ^= self.gf.gf_mult(divisor[j], coef)

                for l in range(len(msg_out) - 1, -1, -1):
                    if msg_out[l] == 0:
                        msg_out.pop()
                    else:
                        break

        dividend_power = len(dividend) - 1
        divisor_power = len(divisor) - 1
        quotient_power = len(quotient) - 1

        if (quotient_power + divisor_power) != dividend_power:
            quotient.insert(0, 0)
        return quotient, msg_out