from re import match

from numpy import float64


SQL_TYPES = {float64 : "REAL",
             int : "INTEGER"}


def ascii_table_records(response, block_size=512):
    """Reads the next line from an ascii table.

    Args:
        response: A http.client.HTTPResponse object.
        block_size: Size in bytes of blocks that will be read.

    Yields:
        Record of the HITRAN database.
    """
    record = ""
    while True:
        block = response.read(block_size).decode("utf-8")
        lines = block.split("\n")
        if lines[-1] == "":
            #If a block ends with a new line character, delete the last
            #element of the list because it will be an empty string.
            del lines[-1]
        for line in lines[:-1]:
            #Return complete lines within a block.
            record += line
            yield record
            record = ""
        if len(block) != block_size:
            #This is the last block.
            yield record + lines[-1]
            break
        elif block.endswith("\n"):
            #No carry-over data between blocks.
            yield lines[-1]
            record = ""
        else:
            #Carry partial last line over to next block.
            record = lines[-1]


def scrub(string):
    """Scrubs user-provided string to prevent database injection.

    Args:
        string: User provided database input.

    Returns:
        The string scrubbed of any trailing spaces, punctuation, or additional text.
    """
    return match(r"([A-Za-z0-9+_-]+)", string.strip()).group(1)
