from mtphandler.model import Well

# regex patterns
WELL_ID_PATTERN = r"[A-H][0-9]{1,2}"


def xy_to_id(x: int, y: int) -> str:
    """Well coordinates to well ID"""
    return f"{chr(y + 65)}{x+1}"


def id_to_xy(well_id: str) -> tuple[int, int]:
    """Convert well ID (e.g., A01 or A1) to well coordinates"""
    row = ord(well_id[0].upper()) - 65
    column = int(well_id[1:]) - 1
    return column, row

    def get_well(self, id: str) -> Well:
        for well in self.plate.wells:
            if well.id.lower() == id.lower():
                return well

        raise ValueError(f"Well {id} not found")


if __name__ == "__main__":
    print(xy_to_id(0, 0))
    print(id_to_xy("A1"))
    print(id_to_xy("H12"))
    print(id_to_xy("C3"))
    print(id_to_xy("A02"))
    print(id_to_xy("A01"))
    print(id_to_xy("A11"))
