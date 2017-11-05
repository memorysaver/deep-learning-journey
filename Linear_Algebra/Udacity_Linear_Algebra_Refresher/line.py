from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class LineException(Exception):
    """Line Exception base class"""


class NoNonZeroElementsLineException(Exception):
    """No nonzero elements found"""


class Line:
    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2
        if not normal_vector:
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        if not constant_term:
            constant_term = Decimal('0')

        self.normal_vector = normal_vector
        self.constant_term = constant_term
        self.set_basepoint()

    def __str__(self):
        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''
            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'
            if not is_initial_term:
                output += ' '
            if abs(coefficient) != 1:
                output += f'{abs(coefficient)}'
            return output

        n = self.normal_vector
        try:
            initial_index = Line.first_nonzero_index(n.coordinates)
            terms = [write_coefficient(n[i], is_initial_term=(i == initial_index))
                     + f'x_{i+1}'
                     for i in range(self.dimension)
                     if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)
        except NoNonZeroElementsLineException:
            output = '0'

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += f' = {constant}'

        return output

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise NoNonZeroElementsLineException()

    def __eq__(self, line):
        if self.normal_vector.is_zero():
            if not line.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - line.constant_term
                return MyDecimal(diff).is_near_zero()
        elif line.normal_vector.is_zero():
            return False

        if not self.is_parallel_to(line):
            return False

        x0 = self.basepoint
        y0 = line.basepoint
        basepoint_difference = x0 - y0
        n = self.normal_vector
        return basepoint_difference.is_orthogonal_to(n)

    def is_parallel_to(self, line):
        return self.normal_vector.is_parallel_to(line.normal_vector)

    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0'] * self.dimension
            initial_index = Line.first_nonzero_index(n.coordinates)
            initial_coefficient = n[initial_index]
            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)
        except NoNonZeroElementsLineException:
            self.basepoint = None

    def intersection_with(self, line):
        try:
            A, B = self.normal_vector.coordinates
            C, D = line.normal_vector.coordinates
            k1 = self.constant_term
            k2 = line.constant_term
            x_numerator = D*k1 - B*k2
            y_numerator = C*k1 + A*k2
            one_over_denom = Decimal('1')/(A*D - B*C)
            return (Vector([x_numerator, y_numerator])
                    .times_scalar(one_over_denom))
        except ZeroDivisionError:
            return self if self == line else None


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


if __name__ == '__main__':
    line1 = Line(normal_vector=Vector(['4.046', '2.836']),
                 constant_term='1.21')
    line2 = Line(normal_vector=Vector(['10.115', '7.09']),
                 constant_term='3.025')
    print(f'intersection 1: {line1.intersection_with(line2)}')
