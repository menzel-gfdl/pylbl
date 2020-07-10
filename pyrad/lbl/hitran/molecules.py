from html.parser import HTMLParser
from urllib.request import urlopen


class _Parser(HTMLParser):
    """Parses HITRAN molecule ids from the HITRAN website.

    Attributes:
        column: Current column in the html table.
        formula: Chemical formula of the current molecule.
        id: HITRAN integer identifier of the current molecule.
        inside_row: Flag telling if the parser is inside a row of a HTML table.
        inside_table: Flag telling if the parser is inside the correct HTML table.
        molecule_ids: Dictionary mapping chemical formulae to their HITRAN ids.
    """
    def __init__(self):
        HTMLParser.__init__(self)
        self.column = 0
        self.formula = ""
        self.id = ""
        self.inside_row = False
        self.inside_table = False
        self.molecule_ids = {}

    def handle_data(self, data):
        """Stores the current HITRAN id or chemical symbol, based on the parser's current state.

        Args:
            data: The data contained in the current HTML tag.
        """
        if self.inside_row:
            if self.column == 1:
                self.id += data.strip()
            elif self.column == 2:
                self.formula += data.strip()

    def handle_endtag(self, tag):
        """Resets flags and stores data based on the parser's current state.

        Args:
           tag: Current HTML tag.
        """
        if tag == "table" and self.inside_table:
            self.inside_table = False
        if tag == "tr" and self.inside_row and self.column > 0:
            self.molecule_ids[self.formula] = int(self.id)
            self.formula = ""
            self.id = ""
            self.inside_row = False

    def handle_starttag(self, tag, attrs):
        """Sets flags based on the parser's current state.

        Args:
            tag: Current HTML tag.
            attrs: List of tuples containing HTML tag attributes.
        """
        if tag == "table" and ("class", "list-table") in attrs:
            self.inside_table = True
        elif tag == "tr" and self.inside_table:
            self.column = 0
            self.inside_row = True
        elif tag == "td" and self.inside_row:
            self.column += 1


def molecules(url):
    """Creates a dictionary mapping molecular chemical formulae to HITRAN ids.

    Args:
        URL to HITRAN molecules table.

    Returns:
        A dictionary mapping molecular chemical formulae to HITRAN ids.
    """
    parser = _Parser()
    parser.feed(urlopen(url).read().decode("utf-8"))
    return parser.molecule_ids
