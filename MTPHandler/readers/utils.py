# regex patterns
WELL_ID_PATTERN = r"[A-H][0-9]{1,2}"


def xy_to_id(x: int, y: int) -> str:
    """Well coordinates to well ID"""
    return f"{chr(y + 65)}{x+1}"


def id_to_xy(well_id: str) -> tuple[int, int]:
    """Well ID to well coordinates"""
    return int(well_id[1:]) - 1, ord(well_id[0].upper()) - 65
