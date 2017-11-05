from decimal import Decimal, getcontext
from math import acos, pi, sqrt

getcontext().prec = 30

DEGREE_PER_RADIAN = 180.0 / pi


class VectorException(Exception):
    """Vector Exception base class"""


class NormalizeVectorException(VectorException):
    """Zero vector cannot be normalized"""


class CompareAngleException(VectorException):
    """Cannot compare angle with zero vector"""


class ComponentParallelError(VectorException):
    """not a unique parallel component"""


class ComponentOrthogonalError(VectorException):
    """not a unique orthogonal component"""


class CrossProductError(VectorException):
    """cross product can only apply to three dimensions"""


class Vector:
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(self.coordinates)
        except ValueError:
            raise ValueError('The coordinates must be non-empty.')
        except TypeError:
            raise TypeError('The coordinates must be an iterable.')

    def __str__(self):
        rounded_result = [round(x, 3) for x in self.coordinates]
        return f'Vector: {rounded_result}'

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def __add__(self, v):
        return Vector([x + y for x, y in zip(self.coordinates, v.coordinates)])

    def __sub__(self, v):
        return Vector([x - y for x, y in zip(self.coordinates, v.coordinates)])

    def times_scalar(self, c):
        return Vector([x * Decimal(c) for x in self.coordinates])

    def magnitude(self):
        return Decimal(sqrt(sum([x**2 for x in self.coordinates])))

    def normalized(self):
        """Compute the unit vector"""
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0') / magnitude)
        except ZeroDivisionError:
            raise NormalizeVectorException('cannot normalize the zero vector')

    def dot(self, v):
        """Compute the inner product"""
        return sum([x * y for x, y in zip(self.coordinates, v.coordinates)])

    def angle_with(self, v, in_degrees=False):
        """Get the angle with two vectors

        The angle in radians of two vectors is the arc cosine of the
        dot product of two vector's unit vector.
        """
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            angle_in_radians = acos(self._clean_cos(u1.dot(u2)))
            if not in_degrees:
                return angle_in_radians
            return angle_in_radians * DEGREE_PER_RADIAN

        except NormalizeVectorException as ex:
            raise CompareAngleException('cannot get angle with zero vector')

    def _clean_cos(self, cos):
        """This method is to clean up precision problem"""
        return min(1, max(-1, cos))

    def is_zero(self, tolerance=1e-10):
        return (self.magnitude() < tolerance)

    def is_parallel_to(self, v):
        return (self.is_zero()
                or v.is_zero()
                or self.angle_with(v) == pi
                or self.angle_with(v) == 0)

    def is_orthogonal_to(self, v, tolerance=1e-10):
        return (abs(self.dot(v)) < tolerance)

    def component_parallel_to(self, basis):
        try:
            normalized_b = basis.normalized()
            return normalized_b.times_scalar(self.dot(normalized_b))
        except NormalizeVectorException:
            raise ComponentParallelError('not a unique parallel component')

    def component_orthogonal_to(self, basis):
        try:
            return (self - self.component_parallel_to(basis))
        except NormalizeVectorException:
            raise ComponentOrthogonalError('not a unique orthogonal component')

    def cross(self, v):
        """Compute cross product of vectors

        Cross product only exists in three dimensions. V x W orthogonal to
        both V and W in 3 dimensions.
        """
        if self.dimention != 3:
            raise CrossProductError('only defines in 3 dimensions.')
        if self.is_parallel_to(v):
            return Vector([0, 0, 0])
        else:
            x = (self.coordinates[1] * v.coordinates[2]
                 - v.coordinates[1] * self.coordinates[2])
            y = -(self.coordinates[0] * v.coordinates[2]
                  - v.coordinates[0] * self.coordinates[2])
            z = (self.coordinates[0] * v.coordinates[1]
                 - v.coordinates[0] * self.coordinates[1])
            return Vector([x, y, z])

    def area_of_triangle_with(self, v):
        return self.area_of_parallelogram_with(v) / Decimal('2.0')

    def area_of_parallelogram_with(self, v):
        return self.cross(v).magnitude()


if __name__ == '__main__':

    # section 2
    my_vector = Vector([1, 2, 3])
    print(my_vector)
    my_vector2 = Vector([1, 2, 3])
    my_vector3 = Vector([-1, 2, 3])
    print(my_vector == my_vector2)
    print(my_vector2 == my_vector3)

    # section 3
    v = Vector([8.218, -9.341]) + Vector([-1.129, 2.111])
    print(v)
    v2 = Vector([7.119, 8.215]) - Vector([-8.223, 0.878])
    print(v2)
    v3 = Vector([1.671, -1.012, -0.318]).times_scalar(7.41)
    print(v3)

    # section 4
    v = Vector([-0.221, 7.437])
    print(f'{v.magnitude()}')
    v = Vector([8.813, -1.331, -6.247])
    print(f'{v.magnitude()}')
    v = Vector([5.581, -2.136])
    print(f'{v.normalized()}')
    v = Vector([1.996, 3.108, -4.554])
    print(f'{v.normalized()}')

    # section 5
    v = Vector([7.887, 4.138])
    w = Vector([-8.802, 6.776])
    print(v.dot(w))
    v = Vector([-5.955, -4.904, -1.874])
    w = Vector([-4.496, -8.755, 7.103])
    print(v.dot(w))
    v = Vector(['3.183', '-7.627'])
    w = Vector(['-2.668', '5.319'])
    print(v.angle_with(w))
    v = Vector(['7.35', '0.221', '5.188'])
    w = Vector(['2.751', '8.259', '3.985'])
    print(v.angle_with(w, True))

    # section 6
    v = Vector(['-7.579', '-7.88'])
    w = Vector(['22.737', '23.64'])
    print(f'{v} is parallel to {w}: {v.is_parallel_to(w)}')
    print(f'{v} is orthogonal to {w}: {v.is_orthogonal_to(w)}')
    v = Vector(['-2.029', '9.97', '4.172'])
    w = Vector(['-9.231', '-6.639', '-7.245'])
    print(f'{v} is parallel to {w}: {v.is_parallel_to(w)}')
    print(f'{v} is orthogonal to {w}: {v.is_orthogonal_to(w)}')
    v = Vector(['-2.328', '-7.284', '-1.214'])
    w = Vector(['-1.821', '1.072', '-2.94'])
    print(f'{v} is parallel to {w}: {v.is_parallel_to(w)}')
    print(f'{v} is orthogonal to {w}: {v.is_orthogonal_to(w)}')
    v = Vector(['2.118', '4.827'])
    w = Vector(['0', '0'])
    print(f'{v} is parallel to {w}: {v.is_parallel_to(w)}')
    print(f'{v} is orthogonal to {w}: {v.is_orthogonal_to(w)}')

    # section 7
    v = Vector(['3.039', '1.879'])
    b = Vector(['0.825', '2.036'])
    parallel_component = v.component_parallel_to(b)
    print(parallel_component)
    v = Vector(['-9.88', '-3.264', '-8.159'])
    b = Vector(['-2.155', '-9.353', '-9.473'])
    orthogonal_component = v.component_orthogonal_to(b)
    print(orthogonal_component)
    v = Vector(['3.009', '-6.172', '3.692', '-2.51'])
    b = Vector(['6.404', '-9.144', '2.759', '8.718'])
    parallel_component = v.component_parallel_to(b)
    orthogonal_component = v.component_orthogonal_to(b)
    print(f'{v} = {parallel_component} + {orthogonal_component}')

    # section 8
    v = Vector(['8.462', '7.893', '-8.187'])
    w = Vector(['6.984', '-5.975', '4.778'])
    cross = v.cross(w)
    print(f'cross product: {cross}')
    v = Vector(['-8.987', '-9.838', '5.031'])
    w = Vector(['-4.268', '-1.861', '-8.866'])
    area_of_parallelogram = v.area_of_parallelogram_with(w)
    print(f'the area of parallelogram: {area_of_parallelogram}')
    v = Vector(['1.5', '9.547', '3.691'])
    w = Vector(['-6.007', '0.124', '5.772'])
    area_of_triangle = v.area_of_triangle_with(w)
    print(f'the area of triangle: {area_of_triangle}')
