"""Periodic table data: symbols, covalent radii, colors, masses."""

from __future__ import annotations

# Index = atomic number; index 0 is the dummy atom "X".
SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
]

_SYMBOL_TO_Z = {s: z for z, s in enumerate(SYMBOLS)}

# Covalent radii in Angstrom (Cordero et al., Dalton Trans. 2008).
_COVALENT_RADII = {
    1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66,
    9: 0.57, 10: 0.58, 11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07,
    16: 1.05, 17: 1.02, 18: 1.06, 19: 2.03, 20: 1.76, 21: 1.70, 22: 1.60,
    23: 1.53, 24: 1.39, 25: 1.39, 26: 1.32, 27: 1.26, 28: 1.24, 29: 1.32,
    30: 1.22, 31: 1.22, 32: 1.20, 33: 1.19, 34: 1.20, 35: 1.20, 36: 1.16,
    37: 2.20, 38: 1.95, 39: 1.90, 40: 1.75, 41: 1.64, 42: 1.54, 43: 1.47,
    44: 1.46, 45: 1.42, 46: 1.39, 47: 1.45, 48: 1.44, 49: 1.42, 50: 1.39,
    51: 1.39, 52: 1.38, 53: 1.39, 54: 1.40, 55: 2.44, 56: 2.15, 57: 2.07,
    58: 2.04, 59: 2.03, 60: 2.01, 61: 1.99, 62: 1.98, 63: 1.98, 64: 1.96,
    65: 1.94, 66: 1.92, 67: 1.92, 68: 1.89, 69: 1.90, 70: 1.87, 71: 1.87,
    72: 1.75, 73: 1.70, 74: 1.62, 75: 1.51, 76: 1.44, 77: 1.41, 78: 1.36,
    79: 1.36, 80: 1.32, 81: 1.45, 82: 1.46, 83: 1.48, 84: 1.40, 85: 1.50,
    86: 1.50, 87: 2.60, 88: 2.21, 89: 2.15, 90: 2.06, 91: 2.00, 92: 1.96,
    93: 1.90, 94: 1.87,
}

# Jmol/CPK element colors as hex RGB.
_CPK_HEX = {
    1: "FFFFFF", 2: "D9FFFF", 3: "CC80FF", 4: "C2FF00", 5: "FFB5B5",
    6: "909090", 7: "3050F8", 8: "FF0D0D", 9: "90E050", 10: "B3E3F5",
    11: "AB5CF2", 12: "8AFF00", 13: "BFA6A6", 14: "F0C8A0", 15: "FF8000",
    16: "FFFF30", 17: "1FF01F", 18: "80D1E3", 19: "8F40D4", 20: "3DFF00",
    21: "E6E6E6", 22: "BFC2C7", 23: "A6A6AB", 24: "8A99C7", 25: "9C7AC7",
    26: "E06633", 27: "F090A0", 28: "50D050", 29: "C88033", 30: "7D80B0",
    31: "C28F8F", 32: "668F8F", 33: "BD80E3", 34: "FFA100", 35: "A62929",
    36: "5CB8D1", 37: "702EB0", 38: "00FF00", 39: "94FFFF", 40: "94E0E0",
    41: "73C2C9", 42: "54B5B5", 43: "3B9E9E", 44: "248F8F", 45: "0A7D8C",
    46: "006985", 47: "C0C0C0", 48: "FFD98F", 49: "A67573", 50: "668080",
    51: "9E63B5", 52: "D47A00", 53: "940094", 54: "429EB0", 55: "57178F",
    56: "00C900", 57: "70D4FF", 58: "FFFFC7", 59: "D9FFC7", 60: "C7FFC7",
    61: "A3FFC7", 62: "8FFFC7", 63: "61FFC7", 64: "45FFC7", 65: "30FFC7",
    66: "1FFFC7", 67: "00FF9C", 68: "00E675", 69: "00D452", 70: "00BF38",
    71: "00AB24", 72: "4DC2FF", 73: "4DA6FF", 74: "2194D6", 75: "267DAB",
    76: "266696", 77: "175487", 78: "D0D0E0", 79: "FFD123", 80: "B8B8D0",
    81: "A6544D", 82: "575961", 83: "9E4FB5", 84: "AB5C00", 85: "754F45",
    86: "428296", 87: "420066", 88: "007D00", 89: "70ABFA", 90: "00BAFF",
    91: "00A1FF", 92: "008FFF",
}

# Standard atomic weights (abridged, common elements).
_MASSES = {
    1: 1.008, 2: 4.003, 3: 6.94, 4: 9.012, 5: 10.81, 6: 12.011, 7: 14.007,
    8: 15.999, 9: 18.998, 10: 20.180, 11: 22.990, 12: 24.305, 13: 26.982,
    14: 28.085, 15: 30.974, 16: 32.06, 17: 35.45, 18: 39.948, 19: 39.098,
    20: 40.078, 21: 44.956, 22: 47.867, 23: 50.942, 24: 51.996, 25: 54.938,
    26: 55.845, 27: 58.933, 28: 58.693, 29: 63.546, 30: 65.38, 31: 69.723,
    32: 72.630, 33: 74.922, 34: 78.971, 35: 79.904, 36: 83.798, 37: 85.468,
    38: 87.62, 39: 88.906, 40: 91.224, 41: 92.906, 42: 95.95, 43: 97.0,
    44: 101.07, 45: 102.906, 46: 106.42, 47: 107.868, 48: 112.414,
    49: 114.818, 50: 118.710, 51: 121.760, 52: 127.60, 53: 126.904,
    54: 131.293, 55: 132.905, 56: 137.327, 57: 138.905, 72: 178.49,
    73: 180.948, 74: 183.84, 75: 186.207, 76: 190.23, 77: 192.217,
    78: 195.084, 79: 196.967, 80: 200.592, 81: 204.38, 82: 207.2,
    83: 208.980, 92: 238.029,
}


def normalize_symbol(value) -> str:
    """Return a canonical element symbol from a symbol string or atomic number."""
    s = str(value).strip()
    if not s:
        return "X"
    if s.lstrip("-").isdigit():
        z = int(s)
        return SYMBOLS[z] if 0 <= z < len(SYMBOLS) else "X"
    letters = "".join(ch for ch in s if ch.isalpha())
    if not letters:
        return "X"
    two = letters[:2].capitalize()
    if len(letters) >= 2 and two in _SYMBOL_TO_Z:
        return two
    one = letters[0].upper()
    return one if one in _SYMBOL_TO_Z else "X"


def symbol_to_z(symbol: str) -> int:
    return _SYMBOL_TO_Z.get(normalize_symbol(symbol), 0)


def covalent_radius(z: int) -> float:
    return _COVALENT_RADII.get(int(z), 1.50)


def cpk_color(z: int) -> tuple[float, float, float]:
    h = _CPK_HEX.get(int(z))
    if h is None:
        return (0.75, 0.75, 0.75)
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def mass(z: int) -> float:
    # Fallback is a rough estimate; fine for display purposes only.
    return _MASSES.get(int(z), 2.5 * int(z))
