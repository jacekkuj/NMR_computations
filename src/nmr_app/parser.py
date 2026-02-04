from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, List


@dataclass(frozen=True)
class NMRParseResult:
    atom_numbers: List[int]
    isotropic_shieldings: List[float]

    def computed_shifts(self, h_ref_ppm: float) -> List[float]:
        # δ(1H) = H_ref − σ_isotropic
        return [h_ref_ppm - s for s in self.isotropic_shieldings]


# --- exact header as requested ---
_EXACT_HEADER = "SCF GIAO Magnetic shielding tensor (ppm):"

# --- proton isotropic line ---
# Handles both:
#   "  74  H    Isotropic =    24.7551   Anisotropy =     9.8348"
#   "  79  H.   Isotropic =    24.7551   Anisotropy =     9.8348"  (if dot appears)
# Also supports scientific notation and Gaussian D-notation:
#   "Isotropic = 2.3456E+01" or "Isotropic = 2.3456D+01"
_H_LINE = re.compile(
    r"^\s*(\d+)\s+H\.?\s+Isotropic\s*=\s*"
    r"([-+]?\d+(?:\.\d+)?(?:[DEe][-+]?\d+)?)\b"
)


def _is_end_of_giao_section(line: str) -> bool:
    s = line.strip()
    if not s:
        return False

    if s.startswith("End of Minotr"):
        return True
    if s.startswith("Population analysis using the SCF density"):
        return True

    if set(s) == {"*"} and len(s) >= 10:
        return True

    if s.endswith(":") and (
        "Population analysis" in s
        or "Mulliken" in s
        or "Dipole moment" in s
        or "SCF Done" in s
    ):
        return True

    return False


def _to_float(num_str: str) -> float:
    # decimal comma + Gaussian D exponent -> python float
    return float(num_str.replace(",", ".").replace("D", "E").replace("d", "E"))


def parse_gaussian_log_for_1h_isotropic(lines: Iterable[str]) -> NMRParseResult:
    """
    Extract ONLY 1H isotropic shieldings from the LAST occurrence of:
        'SCF GIAO Magnetic shielding tensor (ppm):'
    """

    all_lines = [ln.rstrip("\n") for ln in lines]

    header_idx = None
    for i, ln in enumerate(all_lines):
        if ln.strip() == _EXACT_HEADER:
            header_idx = i

    if header_idx is None:
        raise ValueError(
            f"Header not found: '{_EXACT_HEADER}'.\n"
            "Make sure your Gaussian log contains this exact section title."
        )

    atom_numbers: List[int] = []
    isotropic: List[float] = []

    for ln in all_lines[header_idx + 1 :]:
        if _is_end_of_giao_section(ln):
            break

        m = _H_LINE.match(ln)
        if m:
            atom_numbers.append(int(m.group(1)))
            isotropic.append(_to_float(m.group(2)))

    if not isotropic:
        raise ValueError(
            "No 1H isotropic shieldings found under the requested section.\n"
            f"Expected lines like: '<idx>  H  Isotropic = <value>' after '{_EXACT_HEADER}'."
        )

    return NMRParseResult(atom_numbers=atom_numbers, isotropic_shieldings=isotropic)


def parse_text_blob(blob: str) -> NMRParseResult:
    return parse_gaussian_log_for_1h_isotropic(blob.splitlines())
