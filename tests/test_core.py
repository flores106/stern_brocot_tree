import unittest
from fractions import Fraction

from stern_brocot_tree import SternBrocotTree, Path


class TestPath(unittest.TestCase):
    def test_root_path_is_empty(self):
        tree = SternBrocotTree()
        p = tree.path_to(Fraction(1, 1))
        self.assertEqual(p.runs, ())
        self.assertEqual(len(p), 0)
        self.assertEqual(str(p), "")

    def test_path_runs_are_run_length_encoded(self):
        tree = SternBrocotTree()
        # 2/1 is R from root
        self.assertEqual(tree.path_to(Fraction(2, 1)).runs, (("R", 1),))
        # 1/2 is L from root
        self.assertEqual(tree.path_to(Fraction(1, 2)).runs, (("L", 1),))
        # 3/2 is R then L: path R L
        self.assertEqual(tree.path_to(Fraction(3, 2)).runs, (("R", 1), ("L", 1)))
        # 3/1 is R then R: path R R
        self.assertEqual(tree.path_to(Fraction(3, 1)).runs, (("R", 2),))

    def test_path_str_representation(self):
        tree = SternBrocotTree()
        self.assertEqual(str(tree.path_to(Fraction(3, 1))), "R2")
        self.assertEqual(str(tree.path_to(Fraction(3, 2))), "R1L1")

    def test_path_steps_yields_individual_directions(self):
        tree = SternBrocotTree()
        p = tree.path_to(Fraction(3, 2))
        self.assertEqual(list(p.steps()), ["R", "L"])


class TestRoundTrip(unittest.TestCase):
    def test_round_trip_small_fractions(self):
        tree = SternBrocotTree()
        fracs = [
            Fraction(1, 1), Fraction(1, 2), Fraction(2, 1), Fraction(1, 3),
            Fraction(3, 1), Fraction(2, 3), Fraction(3, 2), Fraction(3, 4),
            Fraction(4, 3), Fraction(5, 7), Fraction(7, 5), Fraction(8, 13),
            Fraction(13, 8), Fraction(21, 34), Fraction(99, 100),
        ]
        for f in fracs:
            with self.subTest(frac=f):
                p = tree.path_to(f)
                self.assertEqual(tree.fraction_at(p), f)

    def test_round_trip_large_numerator_denominator(self):
        tree = SternBrocotTree()
        f = Fraction(999983, 1000003)  # both prime
        p = tree.path_to(f)
        self.assertEqual(tree.fraction_at(p), f)

    def test_round_trip_unreduced_input_is_rejected(self):
        # Fraction auto-reduces, so 2/4 becomes 1/2. We test that the
        # library sees the reduced form.
        tree = SternBrocotTree()
        self.assertEqual(tree.path_to(Fraction(2, 4)), tree.path_to(Fraction(1, 2)))


class TestFractionAt(unittest.TestCase):
    def test_empty_path_is_root(self):
        tree = SternBrocotTree()
        self.assertEqual(tree.fraction_at(Path(runs=())), Fraction(1, 1))

    def test_known_fractions(self):
        tree = SternBrocotTree()
        # R -> 2/1
        self.assertEqual(tree.fraction_at(Path(runs=(("R", 1),))), Fraction(2, 1))
        # L -> 1/2
        self.assertEqual(tree.fraction_at(Path(runs=(("L", 1),))), Fraction(1, 2))
        # R R -> 3/1
        self.assertEqual(tree.fraction_at(Path(runs=(("R", 2),))), Fraction(3, 1))
        # L L -> 1/3
        self.assertEqual(tree.fraction_at(Path(runs=(("L", 2),))), Fraction(1, 3))
        # R L -> 3/2
        self.assertEqual(tree.fraction_at(Path(runs=(("R", 1), ("L", 1)))), Fraction(3, 2))
        # L R -> 2/3
        self.assertEqual(tree.fraction_at(Path(runs=(("L", 1), ("R", 1)))), Fraction(2, 3))


class TestNeighbours(unittest.TestCase):
    def test_root_neighbours(self):
        tree = SternBrocotTree()
        left, right = tree.neighbours(Fraction(1, 1))
        self.assertEqual(left, Fraction(0, 1))
        self.assertEqual(right, float("inf"))

    def test_neighbours_of_two_over_one(self):
        tree = SternBrocotTree()
        # 2/1: left neighbour is 1/1, right is infinity (1/0)
        left, right = tree.neighbours(Fraction(2, 1))
        self.assertEqual(left, Fraction(1, 1))
        self.assertEqual(right, float("inf"))

    def test_neighbours_of_three_over_two(self):
        tree = SternBrocotTree()
        # 3/2: parents in the tree are 1/1 (left) and 2/1 (right)
        left, right = tree.neighbours(Fraction(3, 2))
        self.assertEqual(left, Fraction(1, 1))
        self.assertEqual(right, Fraction(2, 1))

    def test_neighbours_satisfy_determinant_property(self):
        # For any fraction a/b with Farey parents c/d and e/f,
        # a*d - b*c = 1 and e*b - f*a = 1 (or equivalently a*f - b*e = -1).
        tree = SternBrocotTree()
        for f in [Fraction(3, 2), Fraction(5, 7), Fraction(8, 13), Fraction(4, 3)]:
            with self.subTest(frac=f):
                left, right = tree.neighbours(f)
                # left < f < right
                self.assertLess(left, f)
                self.assertLess(f, right)
                # determinant: f.numerator * left.denominator - f.denominator * left.numerator == 1
                self.assertEqual(
                    f.numerator * left.denominator - f.denominator * left.numerator,
                    1,
                )
                # right.numerator * f.denominator - right.denominator * f.numerator == 1
                self.assertEqual(
                    right.numerator * f.denominator - right.denominator * f.numerator,
                    1,
                )


class TestDepth(unittest.TestCase):
    def test_root_depth_zero(self):
        tree = SternBrocotTree()
        self.assertEqual(tree.depth(Fraction(1, 1)), 0)

    def test_depth_of_children(self):
        tree = SternBrocotTree()
        self.assertEqual(tree.depth(Fraction(2, 1)), 1)
        self.assertEqual(tree.depth(Fraction(1, 2)), 1)

    def test_depth_of_grandchildren(self):
        tree = SternBrocotTree()
        self.assertEqual(tree.depth(Fraction(3, 1)), 2)
        self.assertEqual(tree.depth(Fraction(1, 3)), 2)
        self.assertEqual(tree.depth(Fraction(3, 2)), 2)
        self.assertEqual(tree.depth(Fraction(2, 3)), 2)


class TestValidation(unittest.TestCase):
    def test_zero_rejected(self):
        tree = SternBrocotTree()
        with self.assertRaises(ValueError):
            tree.path_to(Fraction(0, 1))

    def test_negative_rejected(self):
        tree = SternBrocotTree()
        with self.assertRaises(ValueError):
            tree.path_to(Fraction(-1, 2))

    def test_non_fraction_rejected(self):
        tree = SternBrocotTree()
        with self.assertRaises(TypeError):
            tree.path_to(0.5)  # type: ignore[arg-type]

    def test_neighbours_rejects_zero(self):
        tree = SternBrocotTree()
        with self.assertRaises(ValueError):
            tree.neighbours(Fraction(0, 1))

    def test_fraction_at_rejects_non_path(self):
        tree = SternBrocotTree()
        with self.assertRaises(TypeError):
            tree.fraction_at("R1")  # type: ignore[arg-type]


class TestPathObject(unittest.TestCase):
    def test_path_equality(self):
        p1 = Path(runs=(("R", 2),))
        p2 = Path(runs=(("R", 2),))
        self.assertEqual(p1, p2)

    def test_path_inequality(self):
        p1 = Path(runs=(("R", 2),))
        p2 = Path(runs=(("R", 1), ("L", 1)))
        self.assertNotEqual(p1, p2)

    def test_path_is_hashable(self):
        p = Path(runs=(("R", 2),))
        s = {p, Path(runs=(("R", 2),))}
        self.assertEqual(len(s), 1)


if __name__ == "__main__":
    unittest.main()
