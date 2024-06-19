
_MARKERS = {
    ".": ".",

    ",": ",",

    "o": "o",
    "circle": "o",
    "round": "o",

    "^": "^",
    "triangle": "^",
    "triangle_up": "^",
    "triangleUp": "^",
    "triangleup": "^",

    "v": "v",
    "triangle_down": "v",
    "triangleDown": "v",
    "triangledown": "v",

    "<": "<",
    "triangle_left": "<",
    "triangleLeft": "<",
    "triangleleft": "<",

    ">": ">",
    "triangle_right": ">",
    "triangleRight": ">",
    "triangleright": ">",

    "2": "2",
    "tri_up": "2",
    "triUp": "2",
    "triup": "2",

    "1": "1",
    "tri_down": "1",
    "triDown": "1",
    "tridown": "1",

    "3": "3",
    "tri_left": "3",
    "triLeft": "3",
    "trileft": "3",

    "4": "4",
    "tri_right": "4",
    "triRight": "4",
    "triright": "4",

    "8": "8",
    "octagon": "8",

    "s": "s",
    "square": "s",

    "p": "p",
    "pentagon": "p",

    "P": "P",
    "plus_filled": "P",
    "plusFilled": "P",
    "plusfilled": "P",

    "*": "*",
    "star": "*",

    "h": "h",
    "hexagon": "h",
    "hexagon1": "h",

    "H": "H",
    "hexagon2": "H",

    "+": "+",
    "plus": "+",
    "plus_thin": "+",
    "plusThin": "+",
    "plusthin": "+",

    "x": "x",
    "x_thin": "x",
    "xThin": "x",
    "xthin": "x",

    "X": "X",
    "x_filled": "X",
    "xFilled": "X",
    "xfilled": "X",

    "D": "D",
    "diamond": "D",

    "d": "d",
    "diamond_thin": "d",
    "diamondThin": "d",
    "diamondthin": "d",
    "thin_diamond": "d",
    "thinDiamond": "d",
    "thindiamond": "d",

    "|": "|",
    "v_line": "|",
    "vLine": "|",
    "vline": "|",

    "_": "_",
    "h_line": "_",
    "hLine": "_",
    "hline": "_",

}


def convert_marker_alias(alias: str) -> str:

    try:
        return _MARKERS[alias]
    except KeyError:
        raise ValueError(f"Unknown marker code or alias: {alias}")
