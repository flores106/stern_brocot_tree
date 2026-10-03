from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterator, List, Optional, Tuple


_INFINITY = float("inf")


def _frac_key(frac):
    """Return a numeric key for comparing fractions, including our infinity sentinel."""
    if frac is _INFINITY:
        return _INFINITY
    return frac


@dataclass(frozen=True)
class Path:
    """A finite path from the root of the Stern-Brocot tree to a fraction.

    The path is a sequence of (direction, run_length) pairs, where direction
    is 'L' (left, smaller) or 'R' (right, larger). A run-length encoding is
    used because every fraction's path is of the form R^a L^b R^c ... with
    at least one step in each run, so this is both compact and canonical.
    The empty path represents the root, 1/1.
    """

    runs: Tuple[Tuple[str, int], ...]

    def __len__(self) -> int:
        return sum(length for _, length in self.runs)

    def steps(self) -> Iterator[str]:
        for direction, length in self.runs:
            for _ in range(length):
                yield direction

    def __str__(self) -> str:
        if not self.runs:
            return ""
        return "".join(f"{d}{n}" for d, n in self.runs)


class SternBrocotTree:
    """Ordered rational retrieval via the Stern-Brocot tree.

    The Stern-Brocot tree enumerates every positive rational in lowest terms
    exactly once. Each node is a fraction a/b; its left child is the mediant
    with its left ancestor and its right child is the mediant with its right
    ancestor. This gives a binary search tree over the rationals ordered by
    value, which is the property this library exploits.

    Design decisions:
    - Fractions are stored as ``fractions.Fraction`` so equality and ordering
      are exact, not floating-point.
    - The path to a fraction is computed by binary search using the mediant
      property, which is O(depth) and depth is O(log) in the numerator and
      denominator for typical fractions.
    - We support any positive rational in lowest terms. Zero and negative
      rationals are out of scope: the Stern-Brocot tree is fundamentally a
      structure for positive rationals, and supporting negatives would require
      a design choice (mirror tree? signed wrapper?) that the brief leaves
      open. We pick the single clear interpretation and reject the rest.
    """

    ROOT = Fraction(1, 1)

    def __init__(self) -> None:
        # The tree is conceptual; we navigate by mediant arithmetic, so no
        # nodes are materialised. This instance exists so users have an object
        # to hang configuration on if the library grows.
        pass

    @staticmethod
    def _validate_fraction(frac: Fraction) -> None:
        if not isinstance(frac, Fraction):
            raise TypeError(f"expected Fraction, got {type(frac).__name__}")
        if frac <= 0:
            raise ValueError(
                f"Stern-Brocot tree only contains positive rationals; got {frac}"
            )
        # Fraction auto-reduces, so we don't need to check gcd here.

    def path_to(self, frac: Fraction) -> Path:
        """Return the path from the root (1/1) to ``frac``.

        Raises ``ValueError`` if ``frac`` is not a positive rational, and
        ``TypeError`` if it isn't a ``Fraction``.
        """
        self._validate_fraction(frac)
        if frac == self.ROOT:
            return Path(runs=())

        # We maintain the enclosing interval (lo, hi) whose mediant is the
        # current node. Initially the whole tree sits between 0/1 and 1/0.
        lo: Fraction = Fraction(0, 1)
        hi = _INFINITY
        current: Fraction = self.ROOT

        steps: List[str] = []

        # The depth of any fraction a/b in the Stern-Brocot tree is bounded;
        # in practice it's O(log(a+b)). We cap iterations to guarantee
        # termination even if input is somehow pathological, though for valid
        # positive Fractions the loop always exits via the equality branch.
        max_iterations = 10 * (frac.numerator + frac.denominator) + 100
        for _ in range(max_iterations):
            if current == frac:
                break
            if _frac_key(frac) < _frac_key(current):
                steps.append("L")
                hi = current
                current = Fraction(lo.numerator + current.numerator, lo.denominator + current.denominator)
            else:
                steps.append("R")
                lo = current
                if hi is _INFINITY:
                    current = Fraction(current.numerator + 1, current.denominator)
                else:
                    current = Fraction(current.numerator + hi.numerator, current.denominator + hi.denominator)
        else:
            # Should be unreachable for valid positive Fractions.
            raise RuntimeError("path search did not converge")

        return Path(runs=self._encode_runs(steps))

    @staticmethod
    def _encode_runs(steps: List[str]) -> Tuple[Tuple[str, int], ...]:
        if not steps:
            return ()
        runs: List[Tuple[str, int]] = []
        prev = steps[0]
        count = 1
        for s in steps[1:]:
            if s == prev:
                count += 1
            else:
                runs.append((prev, count))
                prev = s
                count = 1
        runs.append((prev, count))
        return tuple(runs)

    def fraction_at(self, path: Path) -> Fraction:
        """Return the fraction at ``path``.

        Raises ``ValueError`` if the path contains characters other than
        'L' and 'R'.
        """
        if not isinstance(path, Path):
            raise TypeError(f"expected Path, got {type(path).__name__}")

        lo: Fraction = Fraction(0, 1)
        hi = _INFINITY
        current: Fraction = self.ROOT

        for direction in path.steps():
            if direction == "L":
                hi = current
                current = Fraction(lo.numerator + current.numerator, lo.denominator + current.denominator)
            elif direction == "R":
                lo = current
                if hi is _INFINITY:
                    current = Fraction(current.numerator + 1, current.denominator)
                else:
                    current = Fraction(current.numerator + hi.numerator, current.denominator + hi.denominator)
            else:
                raise ValueError(f"invalid direction {direction!r}; expected 'L' or 'R'")

        return current

    def neighbours(self, frac: Fraction) -> Tuple[Fraction, Fraction]:
        """Return the (left, right) neighbours of ``frac`` in the tree.

        The left neighbour is the largest fraction smaller than ``frac`` that
        would appear as an ancestor during a search; the right neighbour is
        the smallest fraction larger than ``frac``. These are the Farey
        parents of ``frac``.
        """
        self._validate_fraction(frac)
        if frac == self.ROOT:
            return (Fraction(0, 1), _INFINITY)

        lo: Fraction = Fraction(0, 1)
        hi = _INFINITY
        current: Fraction = self.ROOT

        max_iterations = 10 * (frac.numerator + frac.denominator) + 100
        for _ in range(max_iterations):
            if current == frac:
                break
            if _frac_key(frac) < _frac_key(current):
                hi = current
                current = Fraction(lo.numerator + current.numerator, lo.denominator + current.denominator)
            else:
                lo = current
                if hi is _INFINITY:
                    current = Fraction(current.numerator + 1, current.denominator)
                else:
                    current = Fraction(current.numerator + hi.numerator, current.denominator + hi.denominator)
        else:
            raise RuntimeError("neighbour search did not converge")

        return (lo, hi)

    def depth(self, frac: Fraction) -> int:
        """Return the depth of ``frac`` in the tree, with the root at depth 0."""
        return len(self.path_to(frac))
